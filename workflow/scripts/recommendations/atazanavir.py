def _recommendation(pt: str) -> str:
    recommendation = ""
    s1 = "ugt1a1"
    if s1 in pt:
        if pt[s1] == "normal_metabolizer":
            recommendation = "No recommendation"
        elif pt[s1] == "intermediate_metabolizer":
            recommendation = "There is no need to avoid prescribing of atazanavir based on UGT1A1 genetic test result. Inform the patient that some patients stop atazanavir because of jaundice (yellow eyes and skin), but that this patient’s genotype makes this unlikely."
        elif pt[s1] == "poor_metabolizer":
            recommendation = "Consider an alternative agent particularly where jaundice would be of concern to the patient."
        else:
            recommendation = "No recommendation"
    elif s1 not in pt: 
        recommendation = "No recommendation"
    return recommendation