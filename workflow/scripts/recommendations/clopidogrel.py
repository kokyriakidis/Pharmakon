def Clopidogrel_recommendation(pt: str) -> str:
    recommendation = ""
    s1 = "cyp2c19"
    if s1 in pt:
        if pt[s1] == "ultrarapid_metabolizer":
            recommendation = "No recommendation"
        if pt[s1] == "normal_metabolizer":
            recommendation = "No recommendation"
        elif pt[s1] == "intermediate_metabolizer":
            recommendation = "Alternative antiplatelet therapy (if no contraindication); e.g., Prasugrel, Ticagrelor."
        elif pt[s1] == "poor_metabolizer":
            recommendation = "Alternative antiplatelet therapy (if no contraindication); e.g., Prasugrel, Ticagrelor."
        else:
            recommendation = "No recommendation"
    elif s1 not in pt: 
        recommendation = "No recommendation"
    return recommendation
