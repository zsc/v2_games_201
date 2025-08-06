#!/usr/bin/env python3
import os
import requests
from concurrent.futures import ThreadPoolExecutor, as_completed
import time

# Corrected URLs for missing PDFs
missing_pdfs = {
    "lec7": [
        ("https://github.com/taichi-dev/games201/releases/download/lec7/lec-7.pdf", "lec7.pdf"),  # Note: lec-7.pdf not lec7.pdf
    ],
    "lec8": [
        ("https://github.com/taichi-dev/games201/releases/download/lec8/lec8_logistics.pdf", "lec8-logistics.pdf"),  # Note: underscore not hyphen
        ("https://github.com/taichi-dev/games201/releases/download/lec8/lec8_mpm.pdf", "lec8-mpm.pdf"),  # Note: underscore not hyphen
    ],
    "lec10": [
        ("https://github.com/taichi-dev/games201/releases/download/lec10/lec-10.pdf", "lec10.pdf"),  # Note: lec-10.pdf not lec10.pdf
    ]
}

def download_pdf(url, filename, folder):
    """Download a single PDF file"""
    filepath = os.path.join(folder, filename)
    
    try:
        print(f"Downloading {filename}...")
        response = requests.get(url, stream=True, timeout=30)
        response.raise_for_status()
        
        with open(filepath, 'wb') as f:
            for chunk in response.iter_content(chunk_size=8192):
                if chunk:
                    f.write(chunk)
        
        print(f"✓ Downloaded {filename} ({os.path.getsize(filepath) / 1024 / 1024:.1f} MB)")
        return True, filename
    except Exception as e:
        print(f"✗ Failed to download {filename}: {str(e)}")
        return False, filename

def main():
    # Create pdfs directory structure
    pdf_dir = "pdfs"
    
    # Flatten all PDFs into a single list
    all_pdfs = []
    for lecture, pdfs in missing_pdfs.items():
        lecture_dir = os.path.join(pdf_dir, lecture)
        if not os.path.exists(lecture_dir):
            os.makedirs(lecture_dir)
        for url, filename in pdfs:
            all_pdfs.append((url, filename, lecture_dir))
    
    print(f"Downloading {len(all_pdfs)} missing PDFs using 4 threads...\n")
    
    start_time = time.time()
    success_count = 0
    failed_files = []
    
    # Use ThreadPoolExecutor with 4 workers
    with ThreadPoolExecutor(max_workers=4) as executor:
        # Submit all download tasks
        future_to_pdf = {
            executor.submit(download_pdf, url, filename, folder): (url, filename) 
            for url, filename, folder in all_pdfs
        }
        
        # Process completed downloads
        for future in as_completed(future_to_pdf):
            success, filename = future.result()
            if success:
                success_count += 1
            else:
                url, _ = future_to_pdf[future]
                failed_files.append((url, filename))
    
    elapsed_time = time.time() - start_time
    
    print(f"\n{'='*50}")
    print(f"Download completed in {elapsed_time:.1f} seconds")
    print(f"Successfully downloaded: {success_count}/{len(all_pdfs)} files")
    
    if failed_files:
        print(f"\nFailed downloads ({len(failed_files)}):")
        for url, filename in failed_files:
            print(f"  - {filename}: {url}")
    
    print(f"{'='*50}")

if __name__ == "__main__":
    main()