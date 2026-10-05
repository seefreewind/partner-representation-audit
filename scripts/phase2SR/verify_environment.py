import os,sys,json,importlib.metadata as md
os.environ['HF_HUB_OFFLINE']='1';os.environ['TRANSFORMERS_OFFLINE']='1'
method=sys.argv[1]
import torch,transformers,numpy,scipy
if method=='sppider':
    from transformers import EsmModel,EsmTokenizer
    import statsmodels,Bio,matplotlib
    assert not any(d.metadata['Name'].lower()=='esm' for d in md.distributions())
else:
    from transformers import T5EncoderModel,AutoTokenizer
    import pytorch_lightning,torchmetrics,xpair
    assert md.version('xpair')=='0.1.0'
    assert md.version('torch')=='2.7.0'
    assert md.version('pytorch-lightning')=='2.5.1'
print(json.dumps({'environment_import_pass':True,'method':method,'python':sys.version,'torch':torch.__version__,'transformers':transformers.__version__,'offline_flags':True}))
