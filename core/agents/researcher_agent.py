import os
import uuid
from typing import List, Dict, Any, Optional
from pydantic_ai import Agent

from utils.parse_args import args
from utils.config import PROXY_CONFIG
from utils.data_handler import DataHandler

for key, value in PROXY_CONFIG.items():
    os.environ[key] = value

class ResearcherAgent():
    def __init__(self):
        self.agent = Agent(
            args.model_name,
            system_prompt=(
                f'''You are a professional chemist with profond knowledge of chemical synthesis. You are provided with a synthesis plan {args.plan_res}. 
                
                You are an expert chemist specializing in retrosynthetic analysis, with extensive experience in designing multi-step synthetic routes for complex organic molecules.

                You are provided with a detailed retrosynthetic analysis plan for a given target molecule.
                
                Based on this plan, you may choose freely between:
                retrieving relevant reactions from a vector or template database, or
                using a template-free reasoning approach to design the route.
                Your first task is to rewrite the provided plan into the following structured format:

                Strategy Index — a unique identifier for each proposed strategy.
                Critical Reaction — indicate whether this reaction is a key (critical) step in the overall synthesis (Yes / No).
                Reaction Name — specify the name or general class of the reaction (e.g., amide coupling, Suzuki coupling, lactam formation, etc.).
                Resulting Precursors — list the immediate precursors or intermediates generated after this disconnection.
                Reaction Conditions / Notes — summarize the typical reagents, catalysts, solvents, temperature, or any special considerations relevant to this reaction.

                Ensure your output is well-structured, concise, and chemically accurate, suitable for guiding downstream automated retrosynthesis planning.

                Each route should be enclosed in a <ROUTE> block and follow this structured format:

                <Route>
                [
                    {
                        'Molecule set': ['<SMILES of current molecule(s)>'],
                        'Rational': '<Explain why this bond is selected for disconnection>',
                        'Product': ['<SMILES of product molecule(s)>'],
                        'Reaction': ['<reaction equation in SMILES, e.g. Product>>Reactant1.Reactant2>'],
                        'Reactants': ['<List of resulting reactant SMILES>'],
                        'Updated molecule set': ['<Current remaining molecules to analyze>']
                    },
                    {
                        'Molecule set': ['<Next molecule(s) to be analyzed>'],
                        'Rational': '<Explain next disconnection step>',
                        'Product': ['<SMILES>'],
                        'Reaction': ['<SMILES reaction>'],
                        'Reactants': ['<SMILES list>'],
                        'Updated molecule set': ['<SMILES list>']
                    },
                    ...
                ]
                </Route>
                '''
            )
        )

        self.data_handler = DataHandler()







