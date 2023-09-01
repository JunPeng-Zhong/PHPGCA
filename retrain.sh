#!/bin/bash
echo retrain
python run_lgcn_bpr_simgcl_batch.py \
		--epochs 2640 \
		--device cuda:0 \
		--lr 0.01 \
		--hidden-channels 128 \
		--cl-rate 0.01 \
		--eps 0.01 \
		--num-layers 3 \
		--batch-size 8192 \
		--train-mode retrain
