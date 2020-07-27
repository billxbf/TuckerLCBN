from util import *
import sys
import warnings
from tqdm import tqdm as tqdm
warnings.filterwarnings('ignore')


T = pd.read_csv('../data/ProcessedData_Full/regonalStaticT4.csv', header=None)
pairList = []
for i in tqdm(range(61)):
    tmp = []
    corlst1 = []
    corlst2 = []
    for j in range(61):
        corlst1 += [np.corrcoef(T.loc[i], T.loc[61+j])[1,0]]
        corlst2 += [np.corrcoef(T.loc[i], T.loc[122+j])[1,0]]
    tmp += [i, np.argmax(corlst1)+61, np.argmax(corlst2)+122]
    pairList += [tmp]
    
pairList = np.array(pairList)
target = np.arange(183).reshape(61,3, order=1)

acc = (((pairList[:,0] == target[:,0]) & (pairList[:,1] == target[:,1])).sum() +
      ((pairList[:,0] == target[:,0]) & (pairList[:,2] == target[:,2])).sum()) /122

target = pd.DataFrame(target, columns=['Base','Pair1','Pair2'])

target.to_csv("../Results/CorrPairing.csv")
print("Accuracy: {}".format(acc))