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
                f"You are a professional chemist with profond knowledge of chemical synthesis. You are provided with a synthesis plan {args.plan}"
            )
        )

        self.data_handler = DataHandler()


