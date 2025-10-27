## 项目结构：

retro_agent
├── README.md
├── api
│   └── main.py
├── core
│   ├── agents
│   │   ├── __init__.py
│   │   ├── checker_agent.py
│   │   ├── planner_agent.py
│   │   └── researcher_agent.py
│   └── coordinator.py
├── data
│   ├── dataset
│   │   ├── blacklist.csv
│   │   ├── origin_dict.csv
│   │   ├── routes
│   │   ├── train_mol_fp_value_step.pt
│   │   └── val_mol_fp_value_step.pt
│   └── storage
│       ├── history
│       └── memory
├── interface
│   ├── __init__.py
│   ├── __pycache__
│   ├── data_handler.py
│   └── researcher_interface.py
├── models
│   ├── __init__.py
│   ├── mlp_retrosyn 
│   ├── rdchiral
│   ├── saved_models
│   ├── template_free_models 如果匹配不到模板就生成
│   └── value_mlp.py 选分子的value function
└── utils 
    ├── __init__.py
    ├── chem_utils.py
    ├── config.py
    └── parse_args.py


