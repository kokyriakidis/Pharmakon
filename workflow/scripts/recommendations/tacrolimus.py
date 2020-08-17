def _recommendation(pt: str) -> str:
    recommendation = ""
    s1 = ""
    if s1 in pt:
        if pt[s1] == "normal_metabolizer":
            recommendation = "Increase starting dose 1.5 to 2 times recommended starting dose. Total starting dose should not exceed 0.3mg/kg/day. Use therapeutic drug monitoring to guide dose adjustments. Further dose adjustments or selection of alternative therapy may be necessary due to other clinical factors (e.g., medication interactions, or hepatic function)."
        elif pt[s1] == "intermediate_metabolizer":
            recommendation = "Increase starting dose 1.5 to 2 times recommended starting dose. Total starting dose should not exceed 0.3mg/kg/day. Use therapeutic drug monitoring to guide dose adjustments. Further dose adjustments or selection of alternative therapy may be necessary due to other clinical factors (e.g., medication interactions, or hepatic function)."
        elif pt[s1] == "poor_metabolizer":
            recommendation = "Initiate therapy with standard recommended dose. Use therapeutic drug monitoring to guide dose adjustments"
        else:
            recommendation = "No recommendation"
    elif s1 not in pt: 
        recommendation = "No recommendation"
    return recommendation