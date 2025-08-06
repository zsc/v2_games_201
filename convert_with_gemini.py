#!/usr/bin/env python3
import os
import subprocess
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
import time

def convert_pdf_with_gemini(pdf_path):
    """Convert a single PDF using Gemini CLI"""
    pdf_file = Path(pdf_path)
    pdf_name = pdf_file.stem
    lecture_dir = pdf_file.parent.name
    
    # Create output directory
    output_dir = Path("gemini_output") / lecture_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    
    output_file = output_dir / f"{pdf_name}.md"
    
    try:
        print(f"Converting {lecture_dir}/{pdf_name}.pdf...")
        
        # Change to the PDF directory and run gemini
        cmd = [
            "gemini",
            "-p", f"verbatim transcribe the text in {pdf_file.name}",
            "-m", "gemini-2.0-flash",
            "-a"
        ]
        
        # Run the command in the PDF's directory
        result = subprocess.run(
            cmd,
            cwd=pdf_file.parent,
            capture_output=True,
            text=True,
            check=True
        )
        
        # Write output to file
        with open(output_file, 'w', encoding='utf-8') as f:
            f.write(result.stdout)
        
        print(f"✓ Converted {lecture_dir}/{pdf_name}.pdf → {output_file}")
        return True, str(output_file)
        
    except subprocess.CalledProcessError as e:
        print(f"✗ Failed to convert {lecture_dir}/{pdf_name}.pdf: {e.stderr}")
        return False, str(pdf_path)
    except Exception as e:
        print(f"✗ Error converting {lecture_dir}/{pdf_name}.pdf: {str(e)}")
        return False, str(pdf_path)

def main():
    # Find all PDFs
    pdf_files = list(Path("pdfs").rglob("*.pdf"))
    
    print(f"Found {len(pdf_files)} PDFs to convert with Gemini\n")
    
    # Create output directory
    Path("gemini_output").mkdir(exist_ok=True)
    
    start_time = time.time()
    success_count = 0
    failed_files = []
    
    # Process PDFs sequentially to avoid rate limits
    # Using ThreadPoolExecutor with max_workers=1 for easier progress tracking
    with ThreadPoolExecutor(max_workers=1) as executor:
        futures = [executor.submit(convert_pdf_with_gemini, pdf_path) for pdf_path in pdf_files]
        
        for future in as_completed(futures):
            success, file_path = future.result()
            if success:
                success_count += 1
            else:
                failed_files.append(file_path)
    
    elapsed_time = time.time() - start_time
    
    print(f"\n{'='*50}")
    print(f"Conversion completed in {elapsed_time:.1f} seconds")
    print(f"Successfully converted: {success_count}/{len(pdf_files)} files")
    
    if failed_files:
        print(f"\nFailed conversions ({len(failed_files)}):")
        for pdf_path in failed_files:
            print(f"  - {pdf_path}")
    
    print(f"{'='*50}")

if __name__ == "__main__":
    main()