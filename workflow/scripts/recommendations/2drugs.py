def _recommendation(pt: str) -> str:
    recommendation = ""
    s1 = ""
    s2 = ""
    if s1 in pt and s2 in pt:
        if pt[s1] == "ultrarapid_metabolizer" and pt[s2] == "ultrarapid_metabolizer":
            recommendation = ""
        elif pt[s1] == "ultrarapid_metabolizer" and pt[s2] == "rapid_metabolizer":
            recommendation = ""
        elif pt[s1] == "ultrarapid_metabolizer" and pt[s2] == "normal_metabolizer":
            recommendation = ""
        elif pt[s1] == "ultrarapid_metabolizer" and pt[s2] == "intermediate_metabolizer":
            recommendation = ""
        elif pt[s1] == "ultrarapid_metabolizer" and pt[s2] == "poor_metabolizer":
            recommendation = ""
        elif pt[s1] == "normal_metabolizer" and pt[s2] == "ultrarapid_metabolizer":
            recommendation = ""
        elif pt[s1] == "normal_metabolizer" and pt[s2] == "rapid_metabolizer":
            recommendation = ""
        elif pt[s1] == "normal_metabolizer" and pt[s2] == "normal_metabolizer":
            recommendation = ""
        elif pt[s1] == "normal_metabolizer" and pt[s2] == "intermediate_metabolizer":
            recommendation = ""
        elif pt[s1] == "normal_metabolizer" and pt[s2] == "poor_metabolizer":
            recommendation = ""
        elif pt[s1] == "intermediate_metabolizer" and pt[s2] == "ultrarapid_metabolizer":
            recommendation = ""
        elif pt[s1] == "intermediate_metabolizer" and pt[s2] == "rapid_metabolizer":
            recommendation = ""
        elif pt[s1] == "intermediate_metabolizer" and pt[s2] == "normal_metabolizer":
            recommendation = ""
        elif pt[s1] == "intermediate_metabolizer" and pt[s2] == "intermediate_metabolizer":
            recommendation = ""
        elif pt[s1] == "intermediate_metabolizer" and pt[s2] == "poor_metabolizer":
            recommendation = ""
        elif pt[s1] == "poor_metabolizer" and pt[s2] == "ultrarapid_metabolizer":
            recommendation = ""
        elif pt[s1] == "poor_metabolizer" and pt[s2] == "rapid_metabolizer":
            recommendation = ""
        elif pt[s1] == "poor_metabolizer" and pt[s2] == "normal_metabolizer":
            recommendation = ""
        elif pt[s1] == "poor_metabolizer" and pt[s2] == "intermediate_metabolizer":
            recommendation = ""
        elif pt[s1] == "poor_metabolizer" and pt[s2] == "poor_metabolizer":
            recommendation = ""
        else:
            recommendation = "No recommendation"
    elif s1 in pt and s2 not in pt:
        if pt[s1] == "ultrarapid_metabolizer":
            recommendation = ""
        if pt[s1] == "rapid_metabolizer":
            recommendation = ""
        elif pt[s1] == "normal_metabolizer":
            recommendation = ""
        elif pt[s1] == "intermediate_metabolizer":
            recommendation = ""
        elif pt[s1] == "poor_metabolizer":
            recommendation = ""
        else:
            recommendation = "No recommendation"
    elif s1 not in pt and s2 in pt:
        if pt[s2] == "ultrarapid_metabolizer":
            recommendation = ""
        if pt[s2] == "rapid_metabolizer":
            recommendation = ""
        if pt[s2] == "normal_metabolizer":
            recommendation = ""
        elif pt[s2] == "intermediate_metabolizer":
            recommendation = ""
        elif pt[s2] == "poor_metabolizer":
            recommendation = ""
        else:
            recommendation = "No recommendation"
    elif s1 not in pt and s2 not in pt: 
        recommendation = "No recommendation"
    return recommendation
