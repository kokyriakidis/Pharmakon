def _recommendation(pt: str) -> str:
    recommendation = ""
    s1 = ""
    if s1 in pt:
        if pt[s1] == "ultrarapid_metabolizer":
            recommendation = ""
        if pt[s1] == "rapid_metabolizer":
            recommendation = ""
        if pt[s1] == "normal_metabolizer":
            recommendation = ""
        elif pt[s1] == "intermediate_metabolizer":
            recommendation = ""
        elif pt[s1] == "poor_metabolizer":
            recommendation = ""
        else:
            recommendation = "No recommendation"
    elif s1 not in pt: 
        recommendation = "No recommendation"
    return recommendation