# coding=utf-8
import os
import sys
import json
import openai
from pathlib import Path
from pydantic_ai import Agent, RunContext
from pydantic import BaseModel


from utils.config import PROXY_CONFIG
from utils.parse_args import args

for key, value in PROXY_CONFIG.items():
    os.environ[key] = value
os.environ["OPENAI_BASE_URL"] = args.base_url

# api_key = os.environ["OPENAI_API_KEY"]
# client = openai.OpenAI(base_url=args.base_url, api_key=api_key, timeout=180)

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


class PlanOutput(BaseModel):
    strategy_index: str
    core_structure_analysis: str
    strategic_bond_disconnection: str
    reaction_type: str
    reaction_name: str
    description: str

class PlannerAgent:
    def __init__(self, checker_res_path=args.checker_res):
        try:
            with open(checker_res_path, 'r', encoding='utf-8') as file:
                checker_res_content = file.read().strip()
        except FileNotFoundError:
            raise ValueError(f"File not found: {checker_res_path}")
        except Exception as e:
            raise ValueError(f"Error reading file: {str(e)}")
        
        self.agent = Agent(  
            args.model_name,
            instructions=(f'''
                You are an expert chemist performing retrosynthetic analysis. Your task is to design a retrosynthetic strategy for a given target molecule represented in SMILES format, and present your reasoning in the structured format defined by PlanOutput.
                          
                You may use the tool 'find_scaffold' to identify and describe the core structure (scaffold) of the target molecule. This refers to the most characteristic or topologically central part of the molecule—such as a fused ring system, heteroaromatic skeleton, or macrocyclic framework. Your goal is to clarify which part of the molecule defines its chemical identity, as this will guide the selection of key disconnection sites.
                Determine the most strategic bond disconnection to initiate retrosynthesis.This should correspond to the key transformation that most effectively simplifies the molecular framework—for example, breaking a major amide, C–C, or C–O bond rather than performing minor functional-group modifications.The chosen disconnection should meaningfully reduce molecular complexity and reveal plausible synthetic precursors.
                Note: Less stable regions of the molecule are typically easier to modify and thus represent likely points of retrosynthetic cleavage.
                          
                For each proposed disconnection, specify:
                    •	Reaction type (e.g., amide coupling, Suzuki cross-coupling, reductive amination, etc.)
                    •	Reaction name (if known)
                    •	Brief chemical rationale — explain why this reaction and disconnection are chemically reasonable.
                          
                Be concise but chemically precise, focusing on the structural logic of the disconnection rather than listing all possible reactions.Your goal is to produce a clear and interpretable retrosynthetic reasoning step.
                          
                Reference Example: For instance, the scaffold of O=C(N[C@@H]1C(=O)N2[C@@H]1SC([C@@H]2C(=O)O)(C)C)COc1ccccc1is C1C(=O)N2C1SC(C2)C The critical step in retrosynthesis is defined as any modification of the scaffold, and the most unstable part of the scaffold is often the most feasible disconnection site. 
                          
                You must strictly follow the output format and return in JSON format.
                '''), output_type=PlanOutput)
        
        @self.agent.tool
        def find_scaffold(ctx:RunContext, smiles: str, generic: bool = False) -> str:
            """
            Tool to identify molecular scaffold.
            ref: https://docs.datamol.io/stable/tutorials/Scaffolds.html
            paper: https://pubs.acs.org/doi/10.1021/acs.jmedchem.5b01746
            计算分子的 Murcko scaffold
            参数：
            smiles : str
                输入分子的 SMILES 字符串
            generic : bool, optional
                是否计算 generic scaffold（将原子类型简化为碳原子），默认 False
            返回：
            scaffold_smiles : str
                scaffold 的 SMILES 表示；如果解析失败则返回 None
            """
            import datamol as dm
            mol = dm.to_mol(smiles)

            if mol is None:
                print(f"Cannot parse {smiles}")
                return None

            try:
                scaffold = dm.to_scaffold_murcko(mol) # Murcko scaffold 
                return dm.to_smiles(scaffold)
            except Exception as e:
                print(f"Error when computing scaffold: {e}")
                return None
    
    def plan_synthesis(self, smiles: str) -> str:
        prompt = smiles
        try:
            result = self.agent.run_sync(prompt)
            return f"Plan for {smiles}: {result.output}"
        except Exception as e:
            return f"{str(e)}"
    

if __name__ == "__main__":
    if not args.input_smiles:
        raise NameError("Please provide an input molecule SMILES string")

    planner_agent = PlannerAgent()
    result = planner_agent.plan_synthesis(args.input_smiles)

    history_path = args.plan_res
    history = []

    if history_path:
        with open(history_path, "r", encoding="utf-8") as f:
            try:
                history = json.load(f)
            except Exception:
                history = []

    existing_entry = next((entry for entry in history if entry["input_smiles"] == args.input_smiles), None)
    if existing_entry:
        existing_entry["result"] = result
    else:
        history.append({
            "input_smiles": args.input_smiles,
            "result": result
        })

    with open(history_path, "w", encoding="utf-8") as f:
        json.dump(history, f, ensure_ascii=False, indent=2)
    print(result)