#!/usr/bin/env python3
import os
import requests
from concurrent.futures import ThreadPoolExecutor, as_completed
from urllib.parse import urlparse
import time

# PDF URLs organized by lecture
pdf_urls = {
    "lec1": [
        ("https://github.com/taichi-dev/games201/releases/download/lec1/lec1-physics.pdf", "lec1-physics.pdf"),
        ("https://github.com/taichi-dev/games201/releases/download/lec1/lec1-taichi.pdf", "lec1-taichi.pdf"),
    ],
    "lec2": [
        ("https://github.com/taichi-dev/games201/releases/download/lec2/lec2-figures.pdf", "lec2-figures.pdf"),
        ("https://github.com/taichi-dev/games201/releases/download/lec2/lec2-math.pdf", "lec2-math.pdf"),
    ],
    "lec3": [
        ("https://github.com/taichi-dev/games201/releases/download/lec3/lec3-logistics.pdf", "lec3-logistics.pdf"),
        ("https://github.com/taichi-dev/games201/releases/download/lec3/lec3-elasticity.pdf", "lec3-elasticity.pdf"),
        ("https://github.com/taichi-dev/games201/releases/download/lec3/lec3-taichi.pdf", "lec3-taichi.pdf"),
    ],
    "lec4": [
        ("https://github.com/taichi-dev/games201/releases/download/lec4/lec4-logistics.pdf", "lec4-logistics.pdf"),
        ("https://github.com/taichi-dev/games201/releases/download/lec4/lec4-eulerian.pdf", "lec4-eulerian.pdf"),
    ],
    "lec5": [
        ("https://github.com/taichi-dev/games201/releases/download/lec5/lec5-fast-poisson.pdf", "lec5-fast-poisson.pdf"),
    ],
    "lec6": [
        ("https://github.com/taichi-dev/games201/releases/download/lec6/lec6-logistics.pdf", "lec6-logistics.pdf"),
        ("https://github.com/taichi-dev/games201/releases/download/lec6/lec-6.pdf", "lec6-fem-topology.pdf"),
    ],
    "lec7": [
        ("https://github.com/taichi-dev/games201/releases/download/lec7/lec7.pdf", "lec7.pdf"),
    ],
    "lec8": [
        ("https://github.com/taichi-dev/games201/releases/download/lec8/lec8-logistics.pdf", "lec8-logistics.pdf"),
        ("https://github.com/taichi-dev/games201/releases/download/lec8/lec8-mpm.pdf", "lec8-mpm.pdf"),
    ],
    "lec9": [
        ("https://github.com/taichi-dev/games201/releases/download/lec9/lec9.pdf", "lec9.pdf"),
    ],
    "lec10": [
        ("https://github.com/taichi-dev/games201/releases/download/lec10/lec10.pdf", "lec10.pdf"),
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
    # Create pdfs directory
    pdf_dir = "pdfs"
    if not os.path.exists(pdf_dir):
        os.makedirs(pdf_dir)
    
    # Flatten all PDFs into a single list
    all_pdfs = []
    for lecture, pdfs in pdf_urls.items():
        lecture_dir = os.path.join(pdf_dir, lecture)
        if not os.path.exists(lecture_dir):
            os.makedirs(lecture_dir)
        for url, filename in pdfs:
            all_pdfs.append((url, filename, lecture_dir))
    
    print(f"Starting download of {len(all_pdfs)} PDFs using 4 threads...\n")
    
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