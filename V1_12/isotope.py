class Isotope:
    def __init__(self, name, mass, decay, abundance, chain):
        self.name = name
        self.Ms = mass
        self.Lam = decay
        self.Abs = abundance
        self.DecayChain = []

