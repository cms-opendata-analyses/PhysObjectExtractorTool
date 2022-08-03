import json

f = open("nevts_recid.txt", "r")
lines = f.readlines()
f.close()
size=len(lines)

variables=[]
for i in range(0, size):
    aux=lines[i].split()
    variables.append([aux[0], int(aux[1]), aux[2] ])

def findRecid(recid):
    result=" "
    for j in range(0,size):
        if(recid==variables[j][0]):
            result= variables[j]    
    return result      

data = {}
data['data'] =  {}
data['ttbar'] = {}
data['single_top_t_chan'] = {}
data['single_atop_t_chan'] = {}
data['single_top_tW'] = {}
data['wjets'] =  {}

#24119
data['data']['nominal'] = {}
data['data']['nominal']['files']=[]
data['data']['nominal']['files'].append({
    'path': findRecid("24119")[2],
    'nevts': findRecid("24119")[1]
})

#24120
data['data']['nominal']['files'].append({
    'path': findRecid("24120")[2],
    'nevts': findRecid("24120")[1]
})
data['data']['nominal']['nevts_total']= findRecid("24119")[1]+findRecid("24120")[1]

#19980
data['ttbar']['nominal']={}
data['ttbar']['nominal']['files']=[]
data['ttbar']['nominal']['files'].append({
    'path': findRecid("19980")[2],
    'nevts': findRecid("19980")[1]
})
data['ttbar']['nominal']['nevts_total']=findRecid("19980")[1]

#19983
data['ttbar']['scaledown']={}
data['ttbar']['scaledown']['files']=[]
data['ttbar']['scaledown']['files'].append({
    'path': findRecid("19983")[2],
    'nevts': findRecid("19983")[1]
})
data['ttbar']['scaledown']['nevts_total']=findRecid("19983")[1]

#19985
data['ttbar']['scaleup']={}
data['ttbar']['scaleup']['files']=[]
data['ttbar']['scaleup']['files'].append({
    'path': findRecid("19985")[2],
    'nevts': findRecid("19985")[1]
})
data['ttbar']['scaleup']['nevts_total']=findRecid("19985")[1]

#19949
data['ttbar']['ME_var']={}
data['ttbar']['ME_var']['files']=[]
data['ttbar']['ME_var']['files'].append({
    'path': findRecid("19949")[2],
    'nevts': findRecid("19949")[1]
})
data['ttbar']['ME_var']['nevts_total']=findRecid("19949")[1]

#19999
data['ttbar']['PS_var']={}
data['ttbar']['PS_var']['files']=[]
data['ttbar']['PS_var']['files'].append({
    'path': findRecid("19999")[2],
    'nevts': findRecid("19999")[1]
})
data['ttbar']['PS_var']['nevts_total']=findRecid("19999")[1]

#19397
data['single_top_t_chan']['nominal'] = {}
data['single_top_t_chan']['nominal']['files']=[]
data['single_top_t_chan']['nominal']['files'].append({
    'path': findRecid("19397")[2],
    'nevts': findRecid("19397")[1]
})
data['single_top_t_chan']['nominal']['nevts_total']=findRecid("19397")[1]

#19407
data['single_atop_t_chan']['nominal'] = {}
data['single_atop_t_chan']['nominal']['files']=[]
data['single_atop_t_chan']['nominal']['files'].append({
    'path': findRecid("19407")[2],
    'nevts': findRecid("19407")[1]
})
data['single_atop_t_chan']['nominal']['nevts_total']=findRecid("19407")[1]

#19419
data['single_top_tW']['nominal'] = {}
data['single_top_tW']['nominal']['files']=[]
data['single_top_tW']['nominal']['files'].append({
    'path': findRecid("19419")[2],
    'nevts': findRecid("19419")[1]
})

#19412
data['single_top_tW']['nominal']['files'].append({
    'path': findRecid("19412")[2],
    'nevts': findRecid("19412")[1]
})
data['single_top_tW']['nominal']['nevts_total']= findRecid("19419")[1]+findRecid("19412")[1]

#20548
data['wjets']['nominal'] = {}
data['wjets']['nominal']['files']=[]
data['wjets']['nominal']['files'].append({
    'path': findRecid("20548")[2],
    'nevts': findRecid("20548")[1]
})
data['wjets']['nominal']['nevts_total']=findRecid("20548")[1]

with open('ntuples.json', 'w') as files:
    json.dump(data, files)
