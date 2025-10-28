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

import pickle
from syntheseus.search.mol_inventory import SmilesListInventory

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

    def __init__(self,inventory_set:set):
        """
        inventory_set: 外部传入的库存集合
        """

        self.data_handler = DataHandler()
        try:
            self.inventory_set = inventory_set
        except:
            FileExistsError("Failed to load inventory set.")


    def check_validity(self) -> ValidityOutput:
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

    def check_buyability(self) -> BuyabilityOutput:

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
                raise ValueError(f"Invalid SMILES string: {smiles}")
            else:
                print
                buyability_map[smiles] = smiles in inventory_set
        return BuyabilityOutput(buyability=buyability_map)
    

if __name__ == "__main__":
    inventory_path = args.mol_inventory_path

    with open(inventory_path, "rb") as f:
        inventory_data = pickle.load(f)
    if isinstance(inventory_data, SmilesListInventory):
        smiles_data = [molecule.smiles for molecule in inventory_data.purchasable_mols()]
        print(f"Inventory loaded from SmilesListInventory with {len(smiles_data)} entries.")
    elif isinstance(inventory_data, list):
        smiles_data = inventory_data
    else:
        raise TypeError("Inventory data format not recognized.")

    checker = CheckerAgent(set(smiles_data))
    res = checker.check_buyability()
    print(res)