class DataHandler:
    def load_data(self, filepath: str) -> dict:
        import json
        with open(filepath, "r") as f:
            data = json.load(f)
        return data
    
    def load_single_input(self, input:str) -> dict:
        return input