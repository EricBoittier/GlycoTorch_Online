def abbreviationToSugarName(abbreviation):
    if abbreviation == "GAD" or abbreviation == "GCD":
        return "ΔUA"
    if abbreviation == "NAG":
        return "GlcNAc"
    if abbreviation == "IDS":
        return "Ido2S"
    if abbreviation == "SGN":
        return "GlcNS6S"
    if abbreviation == "UAP":
        return "ΔUA4S"
    if abbreviation == "SUS":
        return "GlcNS2S6S"
    if abbreviation == "IDY":
        return "Ido2S"
    if abbreviation == "IDY":
        return "Ido2S"
    if abbreviation == "NGA":
        return "GlcNAc"
    if abbreviation == "GNX":
        return "GlcNS3S"
    if abbreviation == "GNS":
        return "GlcNS"
    if abbreviation == "ASG":
        return "GlcNAc4S"
    if abbreviation == "GC4":
        return "GlcUA"
    if abbreviation == "BDP":
        return "GlcUA"
    if abbreviation == "NGK":
        return "GalNAc4S"
    if abbreviation == "NG6":
        return "GalNAc6S"
    if abbreviation == "IDU":
        return "Ido2S"
    if abbreviation == "IXD":
        return "ΔUA2S"
    if abbreviation == "NGY":
        return "GlcNAc6S"
    if abbreviation == "GCU":
        return "GlcUA"
    else:
        return abbreviation

def getSNFGname(abbreviation):
    if abbreviation == "GAD" or abbreviation == "GCD":
        return "GlcA()"
    if abbreviation == "NAG":
        return "GlcNAc()"
    if abbreviation == "IDS":
        return "GlcA( -D \"2S\")"
    if abbreviation == "SGN":
        return "GlcA( -D \"NS6S\")"
    if abbreviation == "UAP":
        return "ΔUA4S()"
    if abbreviation == "SUS":
        return "GlcA( -D \"NS2S6S\")"
    if abbreviation == "IDY":
        return "Ido( -D \"2S\")"
    if abbreviation == "IDY":
        return "Ido( -D \"2S\")"
    if abbreviation == "NGA":
        return "GlcNAc()"
    if abbreviation == "GNX":
        return "GlcN( -D \"NS3S\")"
    if abbreviation == "GNS":
        return "GlcN( -D \"NS\")"
    if abbreviation == "ASG":
        return "GlcNAc( -D \"2S\")"
    if abbreviation == "GC4":
        return "GlcA()"
    if abbreviation == "BDP":
        return "GlcA()"
    if abbreviation == "NGK":
        return "GalNAc( -D \"2S\")"
    if abbreviation == "NG6":
        return "GalNAc( -D \"6S\")"
    if abbreviation == "IDU":
        return "Ido( -D \"2S\")"
    if abbreviation == "IXD":
        return "GlcA( -D \"2S\")"
    else:
        return abbreviation


def catchNameExceptions(name):
        if name == "H1S-(a1-4)-H1S":
            name = "ΔUA2S-(a1-4)-GlcNS6S"

        if name == "H3S-(a1-4)-H3S":
            name = "ΔUA3S-(a1-4)-GlcNS"

        if name == "UCD-(B1-3)-UCD":
            name = "ΔUA-(B1-3)-GalNAc"

        if name == "L42-(B1-3)-L42":
            name = "ΔUA-(B1-3)-GalNAc1S"

        return name
