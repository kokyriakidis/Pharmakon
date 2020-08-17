def Fluorouracil_recommendation(pt: str) -> str:
    recommendation = ""
    s1 = "dpyd"
    if s1 in pt:
        if pt[s1] == "normal_metabolizer":
            recommendation = "Based on genotype, there is no indication to change dose or therapy. Use label-recommended dosage and administration."
        elif pt[s1] == "intermediate_metabolizer":
            recommendation = "Reduce starting dose by 50% followed by titration of dose based on toxicity or therapeutic drug monitoring."
        elif pt[s1] == "poor_metabolizer":
            recommendation = "Avoid use of 5- fluorouracil or 5-fluorouracil prodrug-based regimens. In the event, based on clinical advice, alternative agents are not considered a suitable therapeutic option, 5-fluorouracil should be administered at a strongly reduced dose with early therapeutic drug monitoring. Therapeutic drug monitoring should be done at the earliest time point possible (e.g., minimum time point in steady state) in order to immediately discontinue the infusion if the drug level is too high."
        else:
            recommendation = "No recommendation"
    elif s1 not in pt: 
        recommendation = "No recommendation"
    return recommendation