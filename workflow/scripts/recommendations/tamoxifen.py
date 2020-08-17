def _recommendation(pt: str) -> str:
    recommendation = ""
    s1 = "cyp2d6"
    if s1 in pt:
        if pt[s1] == "ultrarapid_metabolizer":
            recommendation = "No recommendation"
        if pt[s1] == "normal_metabolizer":
            recommendation = "No recommendation"
        elif pt[s1] == "intermediate_metabolizer":
            recommendation = "This patient is predicted to have lower than normal CYP2D6 activity and may be at an increased risk of a poor tamoxifen response due to low endoxifen (active component of tamoxifen) concentrations. Consider hormonal therapy such as an aromatase inhibitor for postmenopausal women or aromatase inhibitor along with ovarian function suppression in premenopausal women. If aromatase inhibitor use is contraindicated, consideration should be given to use a higher but FDA approved tamoxifen dose (40 mg/day). Avoid CYP2D6 inhibitors. Please consult a clinical pharmacistb for more information."
        elif pt[s1] == "poor_metabolizer":
            recommendation = "This patient is predicted to be a CYP2D6 poor metabolizer and may be at an increased risk of a poor tamoxifen response due to low endoxifen (active component of tamoxifen) concentrations. Recommend alternative hormonal therapy such as an aromatase inhibitor for postmenopausal women or aromatase inhibitor along with ovarian function suppression in premenopausal women. Note, higher dose tamoxifen (40 mg/day) increases but does not normalize endoxifen concentrations in CYP2D6 poor metabolizers."
        else:
            recommendation = "No recommendation"
    elif s1 not in pt: 
        recommendation = "No recommendation"
    return recommendation