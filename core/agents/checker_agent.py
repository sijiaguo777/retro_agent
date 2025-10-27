# coding=utf-8
"""
checker agent
"""

import os
import json

from pathlib import Path
from rdkit import Chem, RDLogger
from pydantic import BaseModel
from pydantic_ai import Agent, RunContext
RDLogger.DisableLog('rdApp.*')

from utils.parse_args import args
from utils.config import PROXY_CONFIG
from utils.data_handler import DataHandler

os.environ["OPENAI_BASE_URL"] = args.base_url
for key, value in PROXY_CONFIG.items():
    os.environ[key] = value


def timed(func):
    import time
    from functools import wraps

    @wraps(func)
    def wrapper(*args, **kwargs):
        start = time.time()
        print(f"Running {func.__name__} ...")
        result = func(*args, **kwargs)
        end = time.time()
        print(f"{func.__name__} completed in {end - start:.2f} seconds.\n")
        return result
    return wrapper


class MoleculeValidity(BaseModel):
    valid: bool
    reason: str | None = None

class ValidityOutput(BaseModel):
    validity: dict[str, MoleculeValidity]

class BuyabilityOutput(BaseModel):
    buyability: dict[str, bool]

class CheckerOutput(BaseModel):
    validity: dict[str, dict]
    buyalidity: dict[str, bool]

class CheckerAgent:

    def __init__(self):
        self.agent = Agent(
            args.model_name,
            system_prompt=(
                "You are a professional chemist with profond knowledge of chemical synthesis. Please check the validity and buyability of the molecule in the list. Make sure to use your chemical knowledge make sure that this molecule is valid, safe and commercially available."
            )
        )
        
        @self.agent.tool
        def check_validity(ctx: RunContext) -> ValidityOutput:
            import csv

            filepath = args.mol_file_path
            input_smiles = args.input_smiles
            blacklist_path = args.mol_blacklist_path

            try:
                if input_smiles:
                    data = self.data_handler.load_single_input(input_smiles)
                else:
                    data = self.data_handler.load_data(filepath)
            except Exception as e:
                raise ValueError(f"Failed to load input data: {str(e)}")
            
            validity_map = {}

            if type(data)==str:
                smiles_list = [data]
            else: 
                if "molecule" in data:
                    smiles_list = [data["molecule"]]
                else:
                    smiles_list = list(data.values())
                
            try:
                blacklist = {}
                with open(blacklist_path) as f:
                    reader = csv.DictReader(f)
                    for row in reader:
                        blacklist[row["smiles"]] = row["reason"]
            except FileNotFoundError:
                print(f"Warning: blacklist file not found at {blacklist_path}, skipping blacklist check.")

            for smiles in smiles_list:
                if isinstance(smiles, (list, tuple)):
                    smiles = smiles[0]
                mol = Chem.MolFromSmiles(smiles)
                if mol is None:
                    validity_map[smiles] = MoleculeValidity(valid=False, reason="invalid SMILES")
                elif smiles in blacklist:
                    validity_map[smiles] = MoleculeValidity(valid=False, reason=blacklist[smiles])
                else:
                    validity_map[smiles] = MoleculeValidity(valid=True, reason=None)
            return ValidityOutput(validity=validity_map)
    
        @self.agent.tool
        def check_buyability(ctx: RunContext) -> BuyabilityOutput:

            filepath = args.mol_file_path
            input_smiles = args.input_smiles
            
            try:
                if input_smiles:
                    data = self.data_handler.load_single_input(input_smiles)
                else:
                    data = self.data_handler.load_data(filepath)
            except Exception as e:
                raise ValueError(f"Failed to load input data: {str(e)}")

            buyability_map = {}

            if type(data)==str:
                smiles_list = [data]
            else: 
                if "molecule" in data:
                    smiles_list = [data["molecule"]]
                else:
                    smiles_list = list(data.values())
                    
            inventory_set = self.inventory_set

            for smiles in smiles_list:
                if isinstance(smiles, (list, tuple)):
                    smiles = smiles[0]
                mol = Chem.MolFromSmiles(smiles)
                if mol is None:
                    buyability_map[smiles] = False
                else:
                    buyability_map[smiles] = smiles in inventory_set
            return BuyabilityOutput(buyability=buyability_map)
        
        self.data_handler = DataHandler()
        self.inventory_set = self._load_inventory(args.mol_inventory_path)

    def _load_inventory(self, inventory_path):
        import pickle
        from syntheseus.search.mol_inventory import SmilesListInventory

        try:
            with open(inventory_path, "rb") as f:
                inventory_data = pickle.load(f)
            if isinstance(inventory_data, SmilesListInventory):
                smiles_data = inventory_data.purchasable_mols()
            elif isinstance(inventory_data, list):
                smiles_data = inventory_data
            else:
                raise TypeError(f"Unexpected data type in {inventory_path}: {type(inventory_data)}")
            print(f"Loaded inventory ({len(smiles_data)} entries) from {inventory_path}")
            return set(smiles_data)
        except Exception as e:
            print(f"Warning: failed to load inventory from {inventory_path}: {e}")
            return set()
            
    

    @timed
    def run_validity(self):
        from pydantic_ai.usage import UsageLimits
        try:
            result = self.agent.run_sync(
                '''
                Please check the molecules using the tool 'check_validity'.Output strictly in the required JSON schema.
                If no molecules are valid, output vide JSON.
                ''', output_type=ValidityOutput)
        except UsageLimits as exc:
            limit = self.agent.usage_limits.request_limit
            self.agent.usage_limits.request_limit = limit * 2
            return self.agent.run_sync("Retry the same task now that the limit is higher.",
            output_type=ValidityOutput)
        checker_res_path = args.checker_res
        try:
            with open(checker_res_path, 'w', encoding='utf-8') as file:
                json.dump(result.output, file, indent=4, ensure_ascii=False)
        except Exception as e:
            raise ValueError(f"Failed to write to {checker_res_path}: {str(e)}")

        return result.output
        
    @timed
    def run_buyability(self):
        from pydantic_ai.usage import UsageLimits
        try:
            result = self.agent.run_sync(
                '''
                Please check the buyability of molecules using the tool 'check_buyability'. Output strictly in the required JSON schema. If no molecules are valid, output vide JSON.
                ''', output_type=BuyabilityOutput,usage_limits=UsageLimits(request_limit=200) )
        except UsageLimits as exc:
            limit = self.agent.usage_limits.request_limit
            self.agent.usage_limits.request_limit = limit * 2
            return self.agent.run_sync("Retry the same task now that the limit is higher.",
            output_type=BuyabilityOutput)
        checker_res_path = args.checker_res

        try:
            with open(checker_res_path, 'w', encoding='utf-8') as file:
                json.dump(result.output, file, indent=4, ensure_ascii=False)
        except Exception as e:
            raise ValueError(f"Failed to write to {checker_res_path}: {str(e)}")
        return result.output

    @timed
    def run_checker(self):
        from pydantic_ai.usage import UsageLimits
        result = self.agent.run_sync(
            '''
            You are a professional chemist. Please first check the validity of input molecules using the tool 'check_validity'. Then check the buyability of valid molecules using the tool 'check_buyability'. You must always output JSON strictly in this format:
            {
            "validity": {
                <SMILES>: {"valid": bool, "reason": str}
                }
            }
            {
                "buyalidity": {
                    <SMILES>: bool
                }
            }
            ''',output_type=CheckerOutput,usage_limits=UsageLimits(request_limit=200))

        checker_res_path = args.checker_res
        try:
            with open(checker_res_path, 'w', encoding='utf-8') as file:
                json.dump(result.output, file, indent=4, ensure_ascii=False)
        except Exception as e:
            raise ValueError(f"Failed to write to {checker_res_path}: {str(e)}")
        return result.output

if __name__ == "__main__":
    checker = CheckerAgent()
    res = checker.run_validity()
    # res = checker.run_checker()
    print(res)