from Bio.PDB.Chain import Chain


class ChainNotFound(KeyError):
    def __init__(self, structure_id: str, model_id: str, chain_id: str):
        message = f'Chain {structure_id}/{model_id}/{chain_id} not in structure'
        super().__init__(message)


class NoResiduesFound(ValueError):
    def __init__(self, chain: Chain) -> None:
        chain_id = chain.get_id()
        model_id = chain.get_parent().get_id()
        structure_id = chain.get_parent().get_parent().get_id()
        message = f'No residues matching filter criteria found in {structure_id}/{model_id}/{chain_id}'
        super().__init__(message)
