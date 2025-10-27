import argparse

parser = argparse.ArgumentParser()

# ===================== hyper_param ===================== #
parser.add_argument('--num_strategies', type=int, default=3)
parser.add_argument('--model_name', type=str, default='gpt-5-mini-2025-08-07')
parser.add_argument('--base_url', type=str, default='https://api.xi-ai.cn/v1')

# ===================== data ===================== #
parser.add_argument('--input_smiles', type=str, default='')
parser.add_argument('--mol_file_path', type=str, default='data/storage/history/mols.json')
parser.add_argument('--mol_inventory_path', type=str, default='data/dataset/inventory.pkl') # for buyability checking
parser.add_argument('--mol_blacklist_path', type=str, default='data/dataset/blacklist.csv') # for validity checking
parser.add_argument('--db_path', type=str, default='') 

# ===================== prompt ===================== #
parser.add_argument('--checker_prompt', type=str, default='instructions/checker.txt') 
parser.add_argument('--planner_prompt', type=str, default='instructions/planner.txt')

# ===================== intermediate results ===================== #
parser.add_argument('--checker_res', type=str, default='data/storage/history/checker_res.json') # checker result path
parser.add_argument('--plan_res', type=str, default='data/storage/history/plan.json')

# ==================== logging ==================== #
parser.add_argument('--log_level', default='INFO', 
                    choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'],
                    help='Set the logging level')
parser.add_argument('--console_level', default='INFO',
                    choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'],  
                    help='Set the console output logging level')
parser.add_argument('--silent', action='store_true',
                    help='Silent mode - minimal output')

args = parser.parse_args()