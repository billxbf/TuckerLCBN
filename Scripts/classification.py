from util import *
import argparse
import sys
import warnings
warnings.filterwarnings('ignore')


#%%
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--task_type", type=str, default="individual", help="type of classification labels")
    parser.add_argument("--nfold", type=int, default=3, help="Number of fold in cross validatoin")
    parser.add_argument("--region_path", type=str, default="../data/ProcessedData_Full/regionSingVec_full.csv", help="Regional input data path")
    parser.add_argument("--voxel_path", type=str, default="../data/ProcessedData_Full/voxelSingVec_full.csv", help="Voxel input data path")
    parser.add_argument("--out_path", type=str, default="../Results/ModelsAcc.csv", help="output path")
    opt = parser.parse_args()
    print(opt)
    
    
    if opt.task_type == "individual":
        target = list(np.arange(61)) * 3    
    elif opt.task_type == "farmworker":
        target = [0] * 22 + [1] * 39 + [0] * 22 + [1] * 39 + [0] * 22 + [1] * 39
    elif opt.task_type == "tasks":
        target = [0] * 22 + [0] * 39 + [1] * 22 + [1] * 39 + [2] * 22 + [2] * 39   
    else:
        sys.exit("Task type not an option.")

    df_reg = pd.read_csv(opt.region_path, header = None)
    df_reg['target'] = target
    
    df_vox = pd.read_csv(opt.voxel_path, header = None)
    df_vox['target'] = target
    
    df_reg_10 = df_reg[np.arange(10)]
    df_reg_10['target'] = target
    
    df_vox_10 = df_vox[np.arange(10)]
    df_vox_10['target'] = target
    
    origReg = TrainSeries(df_reg, opt.nfold, False)
    origVox = TrainSeries(df_vox, opt.nfold, False)
    smoteReg = TrainSeries(df_reg, opt.nfold, True)
    smoteVox = TrainSeries(df_vox, opt.nfold, True)
    origReg10 = TrainSeries(df_reg_10, opt.nfold, False)
    origVox10 = TrainSeries(df_vox_10, opt.nfold, False)
    
    
    Model_List = ['Dummy', 'RandomForest','GBDT','KNN','Ridge','SVM']
    AccDF = pd.DataFrame({'Region': origReg, 'Vox': origVox, 
                      'Region_Top10U': origReg10, 'Vox_Top10U':  origVox10, 
                      'Region_SMOTE': smoteReg, 'Vox_SMOTE': smoteVox}, index = Model_List)
            
    AccDF.to_csv(opt.out_path)
    
    
