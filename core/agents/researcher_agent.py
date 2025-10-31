import os
import uuid
import json
from typing import List, Dict, Any, Optional
from pydantic_ai import Agent

from utils.parse_args import args
from utils.config import PROXY_CONFIG
from utils.data_handler import DataHandler

os.environ["OPENAI_BASE_URL"] = args.base_url

os.environ["OPENAI_BASE_URL"] = args.base_url
for key, value in PROXY_CONFIG.items():
    os.environ[key] = value


class ResearcherAgent():
    def __init__(self,plan:str,plan_res:str):
        self.agent = Agent(
            args.model_name,
            system_prompt=(
                
                "You are a professional chemist with profond knowledge of chemical synthesis. " +
                "You are provided with a synthesis plan:\n " + str(plan) +
                '''
                You are an expert chemist specializing in retrosynthetic analysis, with extensive experience in designing multi-step synthetic routes for complex organic molecules.

                "You are provided with a detailed retrosynthetic analysis plan for a given target molecule.\n\n"

                "Based on this plan, you may choose freely between:\n"
                "1. Retrieving relevant reactions from a vector or template database, or\n"
                "2. Using a template-free reasoning approach to design the route.\n\n"

                "Your first task is to rewrite the provided plan into the following structured format:\n\n"

                "Strategy Index — a unique identifier for each proposed strategy.\n"
                "Critical Reaction — indicate whether this reaction is a key (critical) step in the overall synthesis (Yes / No).\n"
                "Reaction Name — specify the name or general class of the reaction (e.g., amide coupling, Suzuki coupling, lactam formation, etc.).\n"
                "Resulting Precursors — list the immediate precursors or intermediates generated after this disconnection.\n"
                "Reaction Conditions / Notes — summarize the typical reagents, catalysts, solvents, temperature, or any special considerations relevant to this reaction.\n\n"

                "Ensure your output is well-structured, concise, and chemically accurate, suitable for guiding downstream automated retrosynthesis planning.\n\n"

                "Each route should be enclosed in a <ROUTE> block and follow this structured format:\n\n"

                "<Route>\n"
                "[\n"
                "    {\n"
                "        'Molecule set': ['<SMILES of current molecule(s)>'],\n"
                "        'Rational': '<Explain why this bond is selected for disconnection>',\n"
                "        'Product': ['<SMILES of product molecule(s)>'],\n"
                "        'Reaction': ['<reaction equation in SMILES, e.g. Product>>Reactant1.Reactant2>'],\n"
                "        'Reactants': ['<List of resulting reactant SMILES>'],\n"
                "        'Updated molecule set': ['<Current remaining molecules to analyze>']\n"
                "    },\n"
                "    {\n"
                "        'Molecule set': ['<Next molecule(s) to be analyzed>'],\n"
                "        'Rational': '<Explain next disconnection step>',\n"
                "        'Product': ['<SMILES>'],\n"
                "        'Reaction': ['<SMILES reaction>'],\n"
                "        'Reactants': ['<SMILES list>'],\n"
                "        'Updated molecule set': ['<SMILES list>']\n"
                "    },\n"
                "    ...\n"
                "]\n"
                "</Route>\n"
                '''
            )
        )

        self.data_handler = DataHandler()


if __name__ == "__main__":
    import json
    plan_path = args.plan_res
    with open(plan_path, "r") as f:
        plan = json.load(f)
        if isinstance(plan, list):
            plan_content = plan[0]["result"]
        elif isinstance(plan, dict):
            plan_content = plan["result"]
        else:
            plan_content = {}

    researcher_agent = ResearcherAgent(plan_content)
    response = researcher_agent.agent.run_sync()
    print("Researcher Agent Response:")
    print(response.output) 
    

