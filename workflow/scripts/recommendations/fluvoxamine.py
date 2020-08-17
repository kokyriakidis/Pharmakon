def _recommendation(pt: str) -> str:
    recommendation = ""
    s1 = "cyp2d6"
    if s1 in pt:
        if pt[s1] == "ultrarapid_metabolizer":
            recommendation = "No recommendation"
        if pt[s1] == "normal_metabolizer":
            recommendation = "No recommendation"
        elif pt[s1] == "intermediate_metabolizer":
            recommendation = "No recommendation"
        elif pt[s1] == "poor_metabolizer":
            recommendation = "This patient is predicted to be a CYP2D6 poor metabolizer and may be at an increased risk of an adverse reaction due to elevated fluvoxamine plasma concentrations. Consider a 25-50% reduction of recommended starting dose and titrate to response or use an alternative drug not metabolized by CYP2D6. Please consult a clinical pharmacist for more information."
        else:
            recommendation = "No recommendation"
    elif s1 not in pt: 
        recommendation = "No recommendation"
    return recommendation