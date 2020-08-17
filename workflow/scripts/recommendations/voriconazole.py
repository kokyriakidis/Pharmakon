def _recommendation(pt: str) -> str:
    recommendation = ""
    s1 = "cyp2c19"
    if s1 in pt:
        if pt[s1] == "ultrarapid_metabolizer":
            recommendation = "Choose an alternative agent that is not dependent on CYP2C19 metabolism as primary therapy in lieu of voriconazole. Such agents include isavuconazole, liposomal amphotericin B, and posaconazole. Further dose adjustments or selection of alternative therapy may be necessary due to other clinical factors, such as drug interactions, hepatic function, renal function, species, site of infection, TDM, and comorbidities."
        if pt[s1] == "rapid_metabolizer":
            recommendation = "Choose an alternative agent that is not dependent on CYP2C19 metabolism as primary therapy in lieu of voriconazole. Such agents include isavuconazole, liposomal amphotericin B, and posaconazole. Further dose adjustments or selection of alternative therapy may be necessary due to other clinical factors, such as drug interactions, hepatic function, renal function, species, site of infection, TDM, and comorbidities."
        if pt[s1] == "normal_metabolizer":
            recommendation = "No recommendation"
        elif pt[s1] == "intermediate_metabolizer":
            recommendation = "No recommendation"
        elif pt[s1] == "poor_metabolizer":
            recommendation = "Choose an alternative agent that is not dependent on CYP2C19 metabolism as primary therapy in lieu of voriconazole. Such agents include isavuconazole, liposomal amphotericin B, and posaconazole.f In the event that voriconazole is considered to be the most appropriate agent, based on clinical advice, for a patient with poor metabolizer genotype, voriconazole should be administered at a preferably lower than standard dosage with careful therapeutic drug monitoring. Further dose adjustments or selection of alternative therapy may be necessary due to other clinical factors, such as drug interactions, hepatic function, renal function, species, site of infection, TDM, and comorbidities."
        else:
            recommendation = "No recommendation"
    elif s1 not in pt: 
        recommendation = "No recommendation"
    return recommendation