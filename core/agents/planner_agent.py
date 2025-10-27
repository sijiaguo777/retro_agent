# coding=utf-8
import os
import sys
import json
from pathlib import Path
from pydantic_ai import Agent
from pydantic import BaseModel

from utils.config import PROXY_CONFIG
from utils.parse_args import args

for key, value in PROXY_CONFIG.items():
    os.environ[key] = value
os.environ["OPENAI_BASE_URL"] = args.base_url


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
                You are an expert chemist performing retrosynthetic analysis.
                Based on the validation results: {checker_res_content}, Your task is to plan a retrosynthetic strategy for a given target molecule represented in SMILES format. You will output your reasoning in the structured format defined by `PlanOutput`. As for core_structure_analysis, you should identify and describe the core structure (or scaffold) of the target molecule. This refers to the most characteristic or topologically central part of the molecule — e.g., a fused ring system, heteroaromatic skeleton, or macrocyclic framework. The goal is to clarify which part of the molecule defines its chemical identity and will guide the identification of key disconnection sites.

                As for strategic_bond_disconnection, you should determine the most strategic bond disconnection* for initiating retrosynthesis.  
                This should correspond to the most important reaction that simplifies the molecule  
                — for instance, breaking a key amide, C–C, or C–O bond that defines the main framework,  
                rather than minor functional group modifications.  
                The chosen disconnection should meaningfully reduce molecular complexity  
                and reveal plausible precursors.

                As for reaction_type / reaction_name / description —  
                Specify the general reaction class (e.g., amide coupling, Suzuki cross-coupling,  
                reductive amination, etc.), give a reaction name if known,  
                and briefly explain why this reaction and disconnection are chemically reasonable.

                Be concise but chemically precise.  
                Focus on the structural logic of the disconnection rather than enumerating all possible reactions.  
                Your goal is to produce a clear, interpretable retrosynthetic reasoning step of planning. Output strictly in the required JSON schema.
                '''), output_type=PlanOutput
        )
        
    @timed
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

    history.append({
        "input_smiles": args.input_smiles,
        "result": result
    })

    with open(history_path, "w", encoding="utf-8") as f:
        json.dump(history, f, ensure_ascii=False, indent=2)
    print(result)