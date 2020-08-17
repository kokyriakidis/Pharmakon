def Codeine_recommendation(pt: str) -> str:
    recommendation = ""
    s1 = "cyp2d6"
    if s1 in pt:
        if pt[s1] == "ultrarapid_metabolizer":
            recommendation = "Avoid Codeine use due to potential for toxicity."
        if pt[s1] == "normal_metabolizer":
            recommendation = "Use label recommended age- or weight-specific dosing."
        elif pt[s1] == "intermediate_metabolizer":
            recommendation = "Use label recommended age- or weight-specific dosing. If no response, consider alternative analgesics such as morphine or a non-opioid."
        elif pt[s1] == "poor_metabolizer":
            recommendation = "Avoid Codeine use due to lack of efficacy."
        else:
            recommendation = "No recommendation"
    elif s1 not in pt: 
        recommendation = "No recommendation"
    return recommendation