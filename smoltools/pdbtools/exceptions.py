class ChainNotFound(KeyError):
    def __init__(self, structure_id: str, model_id: str, chain_id: str):
        message = f'{structure_id}/{model_id}/{chain_id} not in structure'
        super().__init__(message)
