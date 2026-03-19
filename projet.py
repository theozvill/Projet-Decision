import random

# Question 1 : Génération de profils de vote
def generation_profile_approbation(n,m,polarisation):
    profile = []
    vote = [random.randint(0,1) for i in range(m)]
    voteOpp = [1-i for i in vote]

    for _ in range(n):
        r = random.random()

        if r < polarisation/2:
            profile.append(voteOpp)

        elif r < polarisation:
            voteModifie = vote.copy()
            i = random.randint(0,m-1)
            voteModifie[i] = 1 - voteModifie[i]
            profile.append(voteModifie)
        
        else:
            profile.append(vote)

    return profile
    

# Question 2 : Génération de profils de vote avec ordre total
def generation_profile_ordretotal(n,m,polarisation):
    profile = []
    vote = list(range(m))
    random.shuffle(vote)
    voteOpp = vote[::-1]
    
    for _ in range(n):
        r = random.random()
        
        if r < polarisation/2:
            profile.append(voteOpp)

        elif r < polarisation:
            voteModifie = vote.copy()
            i = random.randint(0,m-1)
            j = random.randint(0,m-1)
            voteModifie[i], voteModifie[j] = voteModifie[j], voteModifie[i]
            profile.append(voteModifie)
        
        else:
            profile.append(vote)

