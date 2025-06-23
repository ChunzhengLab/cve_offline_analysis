#!/usr/bin/env python3
"""
Test script for feeddown dispose functionality.
"""

import pandas as pd
import numpy as np
import os
import sys
from pathlib import Path

# Add the current directory to path to import feeddown_dispose
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from feeddown_dispose import calculate_fp, feeddown_correction, process_proton_data, get_lambda_values

def create_test_data():
    """
    Create test data to verify feeddown correction functionality.
    """
    # Create test directory structure
    test_proton_dir = Path("test_data/Proton")
    test_lambda_dir = Path("test_data/Lambda")
    test_proton_dir.mkdir(parents=True, exist_ok=True)
    test_lambda_dir.mkdir(parents=True, exist_ok=True)
    
    # Sample Proton data
    proton_data = {
        'centrality': [5.0, 5.0, 5.0, 15.0, 15.0, 15.0],
        'diff_type': ['Intg', 'Intg', 'Intg', 'Intg', 'Intg', 'Intg'],
        'diff_bin': [0.5, 0.5, 0.5, 0.5, 0.5, 0.5],
        'pair_type': ['SS', 'OS', 'Del', 'SS', 'OS', 'Del'],
        'delta': [-0.001, 0.0001, 0.0011, -0.0012, 0.0002, 0.0014],
        'delta_err': [5e-5, 5e-5, 7e-5, 6e-5, 6e-5, 8e-5],
        'gamma': [2e-5, 3e-5, 1e-5, 1.5e-5, 3.5e-5, 2e-5],
        'gamma_err': [7e-5, 7e-5, 1e-4, 8e-5, 8e-5, 1.1e-4]
    }
    
    # Sample Lambda data
    lambda_data = {
        'centrality': [5.0, 5.0, 5.0, 15.0, 15.0, 15.0],
        'pair_type': ['SS', 'OS', 'Del', 'SS', 'OS', 'Del'],
        'delta': [-0.001038, 0.000698, 0.001737, -0.001306, 0.000752, 0.002058],
        'delta_err': [8.7e-5, 8.7e-5, 1.2e-4, 1.1e-4, 1.1e-4, 1.6e-4],
        'gamma': [-0.000208, 8.6e-5, 0.000294, -0.000687, 9.2e-5, 0.000779],
        'gamma_err': [1.2e-4, 1.2e-4, 1.7e-4, 1.3e-4, 1.3e-4, 1.8e-4]
    }
    
    proton_df = pd.DataFrame(proton_data)
    lambda_df = pd.DataFrame(lambda_data)
    
    # Save test files
    proton_df.to_csv(test_proton_dir / "finalise_test.csv", index=False)
    lambda_df.to_csv(test_lambda_dir / "finalise_test.csv", index=False)
    
    return proton_df, lambda_df

def test_calculate_fp():
    """
    Test the fp calculation.
    """
    print("Testing fp calculation...")
    
    fp, fp_err = calculate_fp()
    
    # Expected values
    expected_fp = (0.8489101154204324 + 0.854837091993265) / 2
    expected_fp_err = np.sqrt(0.12451480395843792**2 + 0.1336355938419346**2) / 2
    
    assert abs(fp - expected_fp) < 1e-10, f"fp mismatch: {fp} vs {expected_fp}"
    assert abs(fp_err - expected_fp_err) < 1e-10, f"fp_err mismatch: {fp_err} vs {expected_fp_err}"
    
    print(f"  ✓ fp = {fp:.6f} ± {fp_err:.6f}")

def test_feeddown_correction():
    """
    Test the feeddown correction function.
    """
    print("\nTesting feeddown correction...")
    
    # Test values
    O_data = 0.001
    sigma_D = 5e-5
    O_lpsec = 0.0005
    sigma_S = 8e-5
    fp = 0.85
    sigma_f = 0.1
    
    O_lp, sigma_Olp = feeddown_correction(O_data, sigma_D, O_lpsec, sigma_S, fp, sigma_f)
    
    # Manual calculation
    expected_O_lp = (O_data - (1 - fp) * O_lpsec) / fp
    
    term1 = sigma_D**2
    term2 = (1 - fp)**2 * sigma_S**2
    term3 = (O_lpsec - expected_O_lp)**2 * sigma_f**2
    expected_sigma = (1 / fp) * np.sqrt(term1 + term2 + term3)
    
    assert abs(O_lp - expected_O_lp) < 1e-10, f"O_lp mismatch: {O_lp} vs {expected_O_lp}"
    assert abs(sigma_Olp - expected_sigma) < 1e-10, f"sigma mismatch: {sigma_Olp} vs {expected_sigma}"
    
    print(f"  ✓ O_lp = {O_lp:.6f} ± {sigma_Olp:.6f}")
    
    # Test NaN handling
    O_lp_nan, sigma_nan = feeddown_correction(np.nan, sigma_D, O_lpsec, sigma_S, fp, sigma_f)
    assert np.isnan(O_lp_nan) and np.isnan(sigma_nan), "Should return NaN for NaN input"
    print("  ✓ NaN handling works correctly")

def test_get_lambda_values():
    """
    Test Lambda value extraction.
    """
    print("\nTesting Lambda value extraction...")
    
    _, lambda_df = create_test_data()
    
    # Test existing values
    delta, delta_err, gamma, gamma_err = get_lambda_values(lambda_df, 5.0, "SS")
    
    expected_row = lambda_df[(lambda_df["centrality"] == 5.0) & (lambda_df["pair_type"] == "SS")].iloc[0]
    
    assert delta == expected_row["delta"], f"Delta mismatch: {delta} vs {expected_row['delta']}"
    assert gamma == expected_row["gamma"], f"Gamma mismatch: {gamma} vs {expected_row['gamma']}"
    
    print(f"  ✓ Found values: delta={delta:.6f}, gamma={gamma:.6f}")
    
    # Test non-existing values
    delta_na, delta_err_na, gamma_na, gamma_err_na = get_lambda_values(lambda_df, 99.0, "SS")
    assert all(np.isnan([delta_na, delta_err_na, gamma_na, gamma_err_na])), "Should return NaN for non-existing data"
    print("  ✓ Non-existing data handling works correctly")

def test_full_processing():
    """
    Test full data processing pipeline.
    """
    print("\nTesting full processing pipeline...")
    
    proton_df, lambda_df = create_test_data()
    fp, fp_err = calculate_fp()
    
    # Process the data
    result_df = process_proton_data(proton_df, lambda_df, fp, fp_err)
    
    if result_df is not None:
        print(f"  ✓ Processing successful. Output shape: {result_df.shape}")
        print(f"  ✓ Columns: {list(result_df.columns)}")
        
        # Check that values have changed (should be different from original)
        original_ss_delta = proton_df[(proton_df['centrality'] == 5.0) & 
                                    (proton_df['pair_type'] == 'SS')]['delta'].iloc[0]
        corrected_ss_delta = result_df[(result_df['centrality'] == 5.0) & 
                                     (result_df['pair_type'] == 'SS')]['delta'].iloc[0]
        
        if abs(corrected_ss_delta - original_ss_delta) > 1e-10:
            print(f"  ✓ Correction applied: {original_ss_delta:.6f} → {corrected_ss_delta:.6f}")
        
        # Check Del recalculation
        intg_group = result_df[(result_df['centrality'] == 5.0) & 
                             (result_df['diff_type'] == 'Intg') & 
                             (result_df['diff_bin'] == 0.5)]
        
        if len(intg_group) == 3:
            ss_delta = intg_group[intg_group['pair_type'] == 'SS']['delta'].iloc[0]
            os_delta = intg_group[intg_group['pair_type'] == 'OS']['delta'].iloc[0]
            del_delta = intg_group[intg_group['pair_type'] == 'Del']['delta'].iloc[0]
            
            expected_del = os_delta - ss_delta
            if abs(del_delta - expected_del) < 1e-10:
                print(f"  ✓ Del recalculation correct: {os_delta:.6f} - {ss_delta:.6f} = {del_delta:.6f}")
        
        # Show sample of processed data
        print("\n  Sample processed data:")
        print(result_df.head(3).to_string(index=False, float_format='%.6f'))
        
    else:
        print("  ✗ Processing failed")

def cleanup_test_data():
    """
    Clean up test data directory.
    """
    import shutil
    if os.path.exists("test_data"):
        shutil.rmtree("test_data")
    print("\n✓ Test data cleaned up")

def main():
    """
    Run all tests.
    """
    print("Feeddown Dispose Test Suite")
    print("=" * 40)
    
    try:
        test_calculate_fp()
        test_feeddown_correction()
        test_get_lambda_values()
        test_full_processing()
        print("\n" + "=" * 40)
        print("All tests passed! ✓")
    except Exception as e:
        print(f"\n✗ Test failed: {e}")
        import traceback
        traceback.print_exc()
    finally:
        cleanup_test_data()

if __name__ == "__main__":
    main()