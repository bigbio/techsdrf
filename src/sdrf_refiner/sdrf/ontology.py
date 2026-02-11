"""
Ontology term mappings for SDRF values.
"""

# Instrument model CV terms (PSI-MS ontology)
INSTRUMENT_TERMS = {
    # Thermo Q Exactive series
    "Q Exactive": "MS:1001911",
    "Q Exactive HF": "MS:1002523",
    "Q Exactive HF-X": "MS:1002877",
    "Q Exactive Plus": "MS:1002634",
    "Exactive": "MS:1000649",
    "Exactive Plus": "MS:1002526",
    # Thermo Orbitrap series
    "Orbitrap Fusion": "MS:1002416",
    "Orbitrap Fusion Lumos": "MS:1002732",
    "Orbitrap Eclipse": "MS:1003029",
    "Orbitrap Exploris 480": "MS:1003028",
    "Orbitrap Exploris 240": "MS:1003096",
    "Orbitrap Exploris 120": "MS:1003095",
    "Orbitrap Astral": "MS:1003378",
    "LTQ Orbitrap": "MS:1000449",
    "LTQ Orbitrap Velos": "MS:1001742",
    "LTQ Orbitrap Elite": "MS:1001910",
    "LTQ Orbitrap XL": "MS:1000556",
    # Thermo LTQ series
    "LTQ": "MS:1000447",
    "LTQ Velos": "MS:1000855",
    "LTQ XL": "MS:1000854",
    # Bruker timsTOF series
    "timsTOF Pro": "MS:1003005",
    "timsTOF Pro 2": "MS:1003230",
    "timsTOF fleX": "MS:1003124",
    "timsTOF SCP": "MS:1003231",
    "timsTOF HT": "MS:1003289",
    "timsTOF Ultra": "MS:1003379",
    # Bruker maXis series
    "maXis": "MS:1000700",
    "maXis II": "MS:1002998",
    "maXis 4G": "MS:1000701",
    # AB SCIEX TripleTOF series
    "TripleTOF 4600": "MS:1002583",
    "TripleTOF 5600": "MS:1000932",
    "TripleTOF 5600+": "MS:1002584",
    "TripleTOF 6600": "MS:1002585",
    "TripleTOF 6600+": "MS:1003047",
    # AB SCIEX ZenoTOF series
    "ZenoTOF 7600": "MS:1003107",
    "ZenoTOF 7600+": "MS:1003377",
    "ZenoTOF 8600": "MS:1003376",
    # Waters
    "Synapt G2": "MS:1001777",
    "Synapt G2-S": "MS:1002088",
    "Synapt G2-Si": "MS:1002726",
    "Xevo G2 QTof": "MS:1001539",
    "Xevo G2-XS QTof": "MS:1002724",
}

# Dissociation method CV terms (PSI-MS ontology)
DISSOCIATION_TERMS = {
    "HCD": "MS:1000422",
    "CID": "MS:1000133",
    "ETD": "MS:1000598",
    "EThcD": "MS:1002631",
    "ETciD": "MS:1002632",
    "UVPD": "MS:1003411",
    "IRMPD": "MS:1000262",
}

# Acquisition method CV terms (PRIDE ontology)
ACQUISITION_TERMS = {
    "Data-dependent acquisition": "PRIDE:0000627",
    "Data-independent acquisition": "PRIDE:0000628",
    "DDA": "PRIDE:0000627",
    "DIA": "PRIDE:0000628",
}


def get_instrument_accession(name: str) -> str:
    """Get CV accession for an instrument name."""
    # Try exact match first
    if name in INSTRUMENT_TERMS:
        return INSTRUMENT_TERMS[name]

    # Try case-insensitive match
    name_lower = name.lower()
    for key, accession in INSTRUMENT_TERMS.items():
        if key.lower() == name_lower:
            return accession

    # Try partial match
    for key, accession in INSTRUMENT_TERMS.items():
        if key.lower() in name_lower or name_lower in key.lower():
            return accession

    return ""


def get_dissociation_accession(name: str) -> str:
    """Get CV accession for a dissociation method."""
    name_upper = name.upper()
    return DISSOCIATION_TERMS.get(name_upper, "")


def get_acquisition_accession(name: str) -> str:
    """Get CV accession for an acquisition method."""
    if name in ACQUISITION_TERMS:
        return ACQUISITION_TERMS[name]

    name_lower = name.lower()
    for key, accession in ACQUISITION_TERMS.items():
        if key.lower() in name_lower:
            return accession

    return ""
