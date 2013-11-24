import os

search_path = r'Z:\data'

LT1_inhomo = []
LT2_inhomo = []

for d, sd, f in os.walk(search_path):
    LT1_inhomo = np.hstack((LT1_inhomo,
        [os.path.join(d,i) for i in f if 'LT1' in i and 'laserscan' in i and 'counts.png' in i and i[38]!='0']))
    LT2_inhomo = np.hstack((LT2_inhomo,
        [os.path.join(d,i) for i in f if 'LT2' in i and 'laserscan' in i and 'counts.png' in i and i[38]!='0']))
