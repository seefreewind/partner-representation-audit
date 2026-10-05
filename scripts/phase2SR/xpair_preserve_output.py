"""Run unchanged official workflow; archive its six-decimal raw export before cleanup."""
from pathlib import Path
import argparse,logging,shutil
from xpair import predict
original=predict.stream_predictions_to_csv
def archive(raw,final,task,threshold):
    shutil.copy2(raw,Path(final).parent/'predictions_six_decimal.csv')
    return original(raw,final,task,threshold)
predict.stream_predictions_to_csv=archive
logging.basicConfig(level=logging.INFO,format='%(asctime)s | %(levelname)s | %(message)s')
params=predict.add_predict_args(argparse.ArgumentParser()).parse_args()
predict.main(params)
