from .constants import CONDA_DEFAULT_ENVIRON
import os

DICT_OF_CONDA_ENVIRONS = {
  "deepLearningToolkit_v0.3": {
    "pythonVersion": "3",
    "dependencies": ["pandas=0.23.4", "scikit-image=0.14.2", "opencv=3.4.2",
                     "tensorflow%(gpuTag)s==1.10.0", "keras=2.2.2",
                     "scikit-learn==0.22"],
    "channels": ["anaconda"],
    "pipPackages": [],
    "defaultInstallOptions": {"gpuTag": ""},  # Tags to be replaced to %(tag)s
    "xmippEnviron": True
  },
  CONDA_DEFAULT_ENVIRON: {
    "pythonVersion": "2.7",
    "dependencies": ["pandas=0.23.4", "scikit-image=0.14.2", "opencv=3.4.2",
                     "tensorflow%(gpuTag)s==1.10.0", "keras=2.2.2",
                     "scikit-learn==0.20"],
    "channels": ["anaconda"],
    "pipPackages": [],
    "defaultInstallOptions": {"gpuTag": ""},  # Tags to be replaced to %(tag)s
    "xmippEnviron": True
  },

  "deepLearningToolkit_v0.01": {
    "pythonVersion": "2.7",
    "dependencies": ["pandas=0.23.4", "scikit-image=0.14.2", "opencv=3.4.2",
                     "tensorflow%(gpuTag)s==1.10.0", "keras=2.1.5",
                     "scikit-learn==0.20"],
    "channels": ["anaconda"],
    "pipPackages": [],
    "defaultInstallOptions": {"gpuTag": ""},  # Tags to be replaced to %(tag)s
    "xmippEnviron": True
  },

  "micrograph_cleaner_em": {
    "pythonVersion": "3.6",
    "dependencies": ["numpy=1.16.4", "micrograph-cleaner-em"],
    "channels": ["rsanchez1369", "anaconda", "conda-forge"],
    "pipPackages": [],
    "defaultInstallOptions": {},
    "xmippEnviron": False
  }

}
