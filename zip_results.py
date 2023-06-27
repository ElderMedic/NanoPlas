import os
import shutil
import zipfile
import pandas as pd

def zip_results(base_dir, sample_config):
    # Group the sample config by Order No.
    grouped_samples = sample_config.groupby('Order No.')
    
    for order_no, samples in grouped_samples:
        # Create the zip file
        zip_filename = os.path.join(base_dir, 'Results', f'{order_no}.zip')
        with zipfile.ZipFile(zip_filename, 'w') as zipf:
            # Iterate over the samples in the current order
            for _, sample in samples.iterrows():
                sample_id = sample['Sample name']
                sample_folder = os.path.join(base_dir, 'Results', sample_id)
                
                # Add each file and folder in the sample folder to the zip file
                for root, _, files in os.walk(sample_folder):
                    for file in files:
                        file_path = os.path.join(root, file)
                        zipf.write(file_path, os.path.relpath(file_path, base_dir))
        
        print(f'Zip file created: {zip_filename}')
