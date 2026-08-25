#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
EVO2 batch mutation effect prediction tool - Advanced version

Optimizations:
- Model is initialized only once
- Smart batching with automatic batch size adjustment
- Checkpoint resume to avoid redundant computation
- GPU memory monitoring
- Error retry mechanism
- Real-time progress saving

Input: mutation_sequences.csv
Output: evo2_predictions.csv (sample_name, Reference_likelihood, Variant_likelihood)
"""

import pandas as pd
import numpy as np
from pathlib import Path
from datetime import datetime
import sys
import gc
import json

# Try importing torch for GPU monitoring
try:
    import torch
    HAS_TORCH = True
except ImportError:
    HAS_TORCH = False

# ==================== Configuration ====================
# Batch size - auto-adjusted based on GPU memory or set manually
INITIAL_BATCH_SIZE = 8      # Initial batch size
MAX_BATCH_SIZE = 32         # Maximum batch size
MIN_BATCH_SIZE = 1          # Minimum batch size

# Memory threshold (GB) - reduce batch size if exceeded
GPU_MEMORY_THRESHOLD = 20

# Retry configuration
MAX_RETRIES = 3             # Maximum number of retries
RETRY_DELAY = 2             # Retry interval (seconds)

# Checkpoint configuration
CHECKPOINT_INTERVAL = 50    # Save checkpoint every N batches
CHECKPOINT_FILE = '.evo2_checkpoint.json'
# ==================================================


def log(msg, level='INFO'):
    """Print a log message with timestamp"""
    timestamp = datetime.now().strftime('%H:%M:%S')
    prefix = {'INFO': '', 'WARN': '⚠️  ', 'ERROR': '❌ ', 'SUCCESS': '✓ '}.get(level, '')
    print(f"[{timestamp}] {prefix}{msg}")


def get_gpu_memory():
    """Get GPU memory usage (GB)"""
    if not HAS_TORCH or not torch.cuda.is_available():
        return None
    try:
        allocated = torch.cuda.memory_allocated() / 1024**3
        reserved = torch.cuda.memory_reserved() / 1024**3
        return {'allocated': allocated, 'reserved': reserved}
    except:
        return None


def log_gpu_memory():
    """Print GPU memory usage"""
    mem = get_gpu_memory()
    if mem:
        log(f"GPU memory: allocated {mem['allocated']:.2f} GB, reserved {mem['reserved']:.2f} GB")


def save_checkpoint(processed_indices, output_data, checkpoint_file=CHECKPOINT_FILE):
    """Save checkpoint"""
    checkpoint = {
        'processed_indices': list(processed_indices) if isinstance(processed_indices, set) else processed_indices,
        'output_data': output_data,
        'timestamp': datetime.now().isoformat()
    }
    with open(checkpoint_file, 'w') as f:
        json.dump(checkpoint, f)


def load_checkpoint(checkpoint_file=CHECKPOINT_FILE):
    """Load checkpoint"""
    if not Path(checkpoint_file).exists():
        return None, []
    try:
        with open(checkpoint_file, 'r') as f:
            checkpoint = json.load(f)
        return set(checkpoint['processed_indices']), checkpoint['output_data']
    except:
        return None, []


def clear_checkpoint(checkpoint_file=CHECKPOINT_FILE):
    """Clear checkpoint"""
    if Path(checkpoint_file).exists():
        Path(checkpoint_file).unlink()


class EVO2Predictor:
    """EVO2 predictor class - encapsulates model and prediction logic"""
    
    def __init__(self, model_name='evo2_7b'):
        """Initialize the EVO2 model"""
        self.model_name = model_name
        self.model = None
        self.batch_size = INITIAL_BATCH_SIZE
        self._initialize_model()
    
    def _initialize_model(self):
        """Load the EVO2 model"""
        log("=" * 70)
        log("Initializing EVO2 model...")
        log(f"Model name: {self.model_name}")
        log("⚠️  Model startup takes about 1 minute, please wait patiently...")
        log("=" * 70)
        
        try:
            from evo2 import Evo2
            self.model = Evo2(self.model_name)
            log("Model initialized successfully!", 'SUCCESS')
            
            # Show model info
            if HAS_TORCH:
                log_gpu_memory()
                
        except ImportError:
            log("Failed to import the evo2 module, please make sure it is installed correctly", 'ERROR')
            raise
        except Exception as e:
            log(f"Model initialization failed: {e}", 'ERROR')
            raise
    
    def _adjust_batch_size(self, success=True):
        """Dynamically adjust batch size based on runtime conditions"""
        if success and self.batch_size < MAX_BATCH_SIZE:
            self.batch_size = min(self.batch_size + 2, MAX_BATCH_SIZE)
        elif not success and self.batch_size > MIN_BATCH_SIZE:
            self.batch_size = max(self.batch_size // 2, MIN_BATCH_SIZE)
    
    def _score_batch(self, sequences, retry=0):
        """
        Score a batch of sequences, with retry mechanism
        
        Args:
            sequences: list of sequences
            retry: current retry count
            
        Returns:
            list: list of scores, returns None on failure
        """
        try:
            scores = self.model.score_sequences(sequences)
            
            # Convert uniformly to a Python list
            if isinstance(scores, torch.Tensor):
                scores = scores.cpu().tolist()
            elif isinstance(scores, np.ndarray):
                scores = scores.tolist()
            elif not isinstance(scores, (list, tuple)):
                scores = [float(scores)]
            
            # Ensure length matches
            if len(scores) != len(sequences):
                log(f"Score count mismatch: expected {len(sequences)}, got {len(scores)}", 'WARN')
                return None
            
            return [float(s) for s in scores]
            
        except RuntimeError as e:
            if "out of memory" in str(e).lower():
                log(f"GPU out of memory, reducing batch size to {self.batch_size // 2}", 'WARN')
                if HAS_TORCH:
                    torch.cuda.empty_cache()
                gc.collect()
                self._adjust_batch_size(success=False)
                
                if retry < MAX_RETRIES:
                    import time
                    time.sleep(RETRY_DELAY)
                    return self._score_batch(sequences, retry + 1)
            raise
            
        except Exception as e:
            log(f"Batch processing failed: {e}", 'WARN')
            if retry < MAX_RETRIES:
                import time
                time.sleep(RETRY_DELAY)
                return self._score_batch(sequences, retry + 1)
            return None
    
    def predict_sequences(self, sequences, desc=""):
        """
        Predict sequences in batch
        
        Args:
            sequences: list of sequences
            desc: description text
            
        Returns:
            list: list of prediction scores
        """
        total = len(sequences)
        all_scores = [None] * total  # Pre-allocate result list
        
        log(f"\nStarting prediction{desc}: {total} sequences")
        log(f"Initial batch size: {self.batch_size}")
        
        idx = 0
        batch_num = 0
        success_count = 0
        
        while idx < total:
            batch_num += 1
            batch_end = min(idx + self.batch_size, total)
            batch = sequences[idx:batch_end]
            
            # Process current batch
            scores = self._score_batch(batch)
            
            if scores is not None:
                all_scores[idx:batch_end] = scores
                success_count += len(batch)
                self._adjust_batch_size(success=True)
            else:
                # Batch failed, process one by one
                log(f"Batch {batch_num} failed, switching to one-by-one processing...", 'WARN')
                for i, seq in enumerate(batch):
                    single_score = self._score_batch([seq])
                    if single_score:
                        all_scores[idx + i] = single_score[0]
                        success_count += 1
            
            idx = batch_end
            
            # Show progress
            if batch_num % 10 == 0 or batch_end >= total:
                progress = batch_end / total * 100
                log(f"  Progress: {batch_end}/{total} ({progress:.1f}%) - "
                    f"success: {success_count}, current batch size: {self.batch_size}")
                
                if batch_num % 50 == 0:
                    log_gpu_memory()
        
        log(f"✓ {desc} prediction complete: {success_count}/{total} succeeded", 'SUCCESS')
        return all_scores
    
    def predict_with_checkpoint(self, df, output_file):
        """
        Prediction with checkpoint, supports resuming from interruption
        
        Args:
            df: input DataFrame
            output_file: output file path
            
        Returns:
            DataFrame: results
        """
        total = len(df)
        
        # Try loading checkpoint
        processed_indices, saved_data = load_checkpoint()
        
        if processed_indices:
            log(f"Checkpoint found, {len(processed_indices)}/{total} already processed, continuing...", 'WARN')
        else:
            processed_indices = set()
            saved_data = []
            log(f"Starting new task: {total} sequences")
        
        # Prepare data
        sample_names = df['sample_name'].tolist()
        wt_sequences = df['wt_sequence'].tolist()
        mut_sequences = df['mut_sequence'].tolist()
        
        ref_scores = [None] * total
        var_scores = [None] * total
        
        # Restore saved data
        for item in saved_data:
            idx = item['index']
            ref_scores[idx] = item['ref_score']
            var_scores[idx] = item['var_score']
        
        # Predict wild-type sequences
        log("\n" + "=" * 70)
        log("Stage 1: Predicting wild-type sequences")
        log("=" * 70)
        
        for i in range(total):
            if i in processed_indices:
                continue
            
            score = self._score_batch([wt_sequences[i]])
            if score:
                ref_scores[i] = score[0]
            
            # Periodically save checkpoint
            if (i + 1) % CHECKPOINT_INTERVAL == 0:
                self._save_progress(i, sample_names, ref_scores, var_scores, output_file)
        
        # Predict mutant sequences
        log("\n" + "=" * 70)
        log("Stage 2: Predicting mutant sequences")
        log("=" * 70)
        
        processed_indices = set()  # Reset for stage 2
        
        for i in range(total):
            if i in processed_indices:
                continue
            
            score = self._score_batch([mut_sequences[i]])
            if score:
                var_scores[i] = score[0]
            
            if (i + 1) % CHECKPOINT_INTERVAL == 0:
                self._save_progress(i, sample_names, ref_scores, var_scores, output_file)
        
        # Build results
        result_df = pd.DataFrame({
            'sample_name': sample_names,
            'Reference_likelihood': ref_scores,
            'Variant_likelihood': var_scores
        })
        
        # Compute difference
        result_df['Likelihood_diff'] = (
            result_df['Variant_likelihood'] - result_df['Reference_likelihood']
        )
        
        # Save final results
        result_df.to_csv(output_file, index=False, encoding='utf-8-sig')
        clear_checkpoint()
        
        return result_df
    
    def _save_progress(self, current_idx, sample_names, ref_scores, var_scores, output_file):
        """Save intermediate progress"""
        output_data = []
        for i in range(current_idx + 1):
            if ref_scores[i] is not None or var_scores[i] is not None:
                # Convert to native Python types to ensure JSON serializability
                ref_val = float(ref_scores[i]) if ref_scores[i] is not None else None
                var_val = float(var_scores[i]) if var_scores[i] is not None else None
                output_data.append({
                    'index': i,
                    'sample_name': sample_names[i],
                    'ref_score': ref_val,
                    'var_score': var_val
                })
        
        save_checkpoint(list(range(current_idx + 1)), output_data)
        log(f"  Checkpoint saved: {current_idx + 1} records")


def main():
    """Main function"""
    # Configuration
    INPUT_FILE = 'mutation_sequences.csv'
    OUTPUT_FILE = 'evo2_predictions.csv'
    MODEL_NAME = 'evo2_7b'
    
    log("=" * 70)
    log("EVO2 batch mutation effect prediction - Advanced version")
    log("=" * 70)
    log(f"Start time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    log(f"PyTorch available: {HAS_TORCH}")
    if HAS_TORCH:
        log(f"CUDA available: {torch.cuda.is_available()}")
    
    # Check input file
    if not Path(INPUT_FILE).exists():
        log(f"Input file does not exist: {INPUT_FILE}", 'ERROR')
        sys.exit(1)
    
    # Read data
    log(f"\nReading input file: {INPUT_FILE}")
    df = pd.read_csv(INPUT_FILE)
    log(f"✓ Total {len(df)} sequences")
    
    # Check column names
    required_cols = ['sample_name', 'wt_sequence', 'mut_sequence']
    for col in required_cols:
        if col not in df.columns:
            log(f"Missing column: {col}", 'ERROR')
            log(f"Available columns: {list(df.columns)}")
            sys.exit(1)
    
    # Initialize predictor
    predictor = EVO2Predictor(MODEL_NAME)
    
    # Run prediction
    result_df = predictor.predict_with_checkpoint(df, OUTPUT_FILE)
    
    # Output statistics
    log("\n" + "=" * 70)
    log("Prediction complete!", 'SUCCESS')
    log("=" * 70)
    log(f"Output file: {OUTPUT_FILE}")
    log(f"Total: {len(result_df)} records")
    
    success_ref = result_df['Reference_likelihood'].notna().sum()
    success_var = result_df['Variant_likelihood'].notna().sum()
    log(f"Reference success: {success_ref} records ({success_ref/len(result_df)*100:.1f}%)")
    log(f"Variant success: {success_var} records ({success_var/len(result_df)*100:.1f}%)")
    
    # Statistics
    log(f"\nStatistics:")
    log(f"  Reference likelihood mean: {result_df['Reference_likelihood'].mean():.4f}")
    log(f"  Variant likelihood mean: {result_df['Variant_likelihood'].mean():.4f}")
    log(f"  Likelihood diff mean: {result_df['Likelihood_diff'].mean():.4f}")
    
    # Preview
    log(f"\nFirst 5 rows preview:")
    print(result_df.head().to_string(index=False))
    
    log(f"\nEnd time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    log("=" * 70)


if __name__ == '__main__':
    main()
