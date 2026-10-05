"""Official unchanged workflow with export preservation and timing only."""
from pathlib import Path
import argparse,logging,shutil,time,json,resource
from xpair import predict
original=predict.stream_predictions_to_csv
forward=predict.XPairModel.predict_step
pairtimes=[]
def timed(self,batch,batch_idx,*args,**kwargs):
    start=time.monotonic();result=forward(self,batch,batch_idx,*args,**kwargs)
    import torch
    if torch.backends.mps.is_available():torch.mps.synchronize()
    elapsed=time.monotonic()-start
    _,data=self._unpack_task_batch(batch)
    for a,b in zip(data['id1'],data['id2']):pairtimes.append(dict(id1=a,id2=b,forward_wall_seconds=elapsed/len(data['id1'])))
    return result
predict.XPairModel.predict_step=timed
def archive(raw,final,task,threshold):
    shutil.copy2(raw,Path(final).parent/'predictions_six_decimal.csv');return original(raw,final,task,threshold)
predict.stream_predictions_to_csv=archive
logging.basicConfig(level=logging.INFO,format='%(asctime)s | %(levelname)s | %(message)s')
params=predict.add_predict_args(argparse.ArgumentParser()).parse_args();start=time.monotonic();cpu=time.process_time();predict.main(params)
(Path(params.predict_output_path).parent/'resource_measurements.json').write_text(json.dumps(dict(pair_timings=pairtimes,total_wall_seconds=time.monotonic()-start,CPU_seconds=time.process_time()-cpu,peak_RAM_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,MPS_allocated_GPU_busy_hours=None,MPS_peak_VRAM=None),indent=2)+'\n')
