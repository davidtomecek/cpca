import os
import glob
import numpy as np
import nibabel as nib

def compute_voxel_fd(func_path, rp_path, output_dir=None):
    """
    Computes 4D voxel-specific Framewise Displacement (FD_vox) and 
    1D whole-brain mean spatial FD according to Yan et al. (2013).
    
    Parameters:
    -----------
    func_path : str
        Path to the 4D functional BOLD NIFTI file (pre-realignment/slice timing).
    rp_path : str
        Path to the SPM-style realignment parameters 'rp_*.txt' (6 columns).
    output_dir : str, optional
        Where to save outputs. Defaults to the directory containing func_path.
    """
    print(f"--- Processing: {os.path.basename(func_path)} ---")
    
    # 1. Load the functional NIFTI image and header affine
    img = nib.load(func_path)
    M_0 = img.affine  # Voxel-to-world mapping matrix (4x4)
    shape = img.shape[:3]
    X_dim, Y_dim, Z_dim = shape
    n_frames = img.shape[3]
    
    # 2. Load SPM realignment parameters (3 translations in mm, 3 rotations in rad)
    # SPM order: X, Y, Z (translations), pitch(alpha), roll(beta), yaw(gamma)
    motion_params = np.loadtxt(rp_path)
    if motion_params.shape[0] != n_frames:
        raise ValueError(f"Mismatch: RP file has {motion_params.shape[0]} frames but NIFTI has {n_frames}.")
        
    if output_dir is None:
        output_dir = os.path.dirname(func_path)
    os.makedirs(output_dir, exist_ok=True)

    # 3. Create a simplistic whole-brain mask based on non-zero background voxels
    # (Replace this with a path to a real skull-stripped mask if available)
    func_data = img.get_fdata()
    brain_mask = np.mean(func_data, axis=-1) > 0
    
    # 4. Construct a 3D coordinate grid of voxel indexes (i, j, k)
    i, j, k = np.meshgrid(np.arange(X_dim), np.arange(Y_dim), np.arange(Z_dim), indexing='ij')
    ones = np.ones_like(i)
    # Flatten to a (4, N_voxels) array for dot-product calculation
    voxel_indices = np.stack([i, j, k, ones], axis=0).reshape(4, -1)
    
    # 5. Transform index grid to baseline real-world millimeter coordinates (C_xyz,0)
    C_xyz_0 = M_0 @ voxel_indices  # Shape: (4, N_voxels)
    
    # Pre-allocate the massive 4D array for Voxel-wise FD map [X x Y x Z x Time]
    fd_vox_4d = np.zeros((X_dim, Y_dim, Z_dim, n_frames), dtype=np.float32)
    
    # Track coordinates over time
    # To compute frame-by-frame changes, we track the coordinates of frame (t) and frame (t-1)
    prev_C_xyz_t = np.copy(C_xyz_0) 
    
    print("Calculating spatial displacements frame by frame...")
    for t in range(n_frames):
        # Extract motion metrics for current frame
        tx, ty, tz, alpha, beta, gamma = motion_params[t, :]
        
        # Build standard SPM rigid-body transformation matrix (T_t)
        # 3D Rotation matrices
        R_x = np.array([[1, 0, 0], [0, np.cos(alpha), -np.sin(alpha)], [0, np.sin(alpha), np.cos(alpha)]])
        R_y = np.array([[np.cos(beta), 0, np.sin(beta)], [0, 1, 0], [-np.sin(beta), 0, np.cos(beta)]])
        R_z = np.array([[np.cos(gamma), -np.sin(gamma), 0], [np.sin(gamma), np.cos(gamma), 0], [0, 0, 1]])
        R = R_z @ R_y @ R_x  # Order of rotations applied by SPM
        
        T_t = np.eye(4)
        T_t[:3, :3] = R
        T_t[:3, 3] = [tx, ty, tz]
        
        # Invert matrix to shift baseline world points back to frame 't' orientation
        T_t_inv = np.linalg.inv(T_t)
        
        # Displaced millimeter coordinates at time t
        C_xyz_t = T_t_inv @ C_xyz_0  # Shape: (4, N_voxels)
        
        if t == 0:
            # The first frame has no backward frame change reference
            fd_vox_frame = np.zeros(C_xyz_0.shape[1])
        else:
            # Calculate absolute Euclidean spatial step distance from (t-1) to (t)
            delta_xyz = C_xyz_t[:3, :] - prev_C_xyz_t[:3, :]
            fd_vox_frame = np.sqrt(np.sum(delta_xyz**2, axis=0))
            
        # Reshape flat voxel row back into standard 3D space matrix coordinates
        fd_vox_4d[:, :, :, t] = fd_vox_frame.reshape(X_dim, Y_dim, Z_dim)
        
        # Cache current frame coordinates as reference for the next cycle loop
        prev_C_xyz_t = C_xyz_t

    # 6. Save the 4D Voxel-Based FD Map to Disk
    sub_id = os.path.basename(rp_path).replace('rp_', '').replace('.txt', '')
    fd_nifti_path = os.path.join(output_dir, f"FDvox_{sub_id}.nii.gz")
    fd_img = nib.Nifti1Image(fd_vox_4d, M_0, img.header)
    nib.save(fd_img, fd_nifti_path)
    print(f"Saved 4D Map: {fd_nifti_path}")
    
    # 7. Calculate and Save Whole-Brain Mean Spatial FD (1D timecourse)
    # Mask out the background to average only inside brain matter space
    mean_sp_fd = np.zeros(n_frames)
    for t in range(n_frames):
        frame_data = fd_vox_4d[:, :, :, t]
        mean_sp_fd[t] = np.mean(frame_data[brain_mask])
        
    mean_fd_path = os.path.join(output_dir, f"mean_sp_FDvox_{sub_id}.txt")
    np.savetxt(mean_fd_path, mean_sp_fd, fmt='%f')
    print(f"Saved 1D Timecourse: {mean_fd_path}\n")
    
    return mean_sp_fd

# --- Pipeline Batch Automation Runner ---
if __name__ == "__main__":
    # Point this to your DPARSF-style dataset layout
    base_dir = "/path/to/your/Study_Working_Dir"
    
    # Loop over detected subjects automatically
    subject_folders = glob.glob(os.path.join(base_dir, "FunImg", "Sub_*"))
    
    for sub_folder in subject_folders:
        sub_id = os.path.basename(sub_folder)
        
        # Find files matching criteria
        func_files = glob.glob(os.path.join(sub_folder, "*.nii*"))
        rp_files = glob.glob(os.path.join(base_dir, "RealignParameter", sub_id, "rp_*.txt"))
        
        if func_files and rp_files:
            # Run the calculation script block
            compute_voxel_fd(
                func_path=func_files[0], 
                rp_path=rp_files[0], 
                output_dir=os.path.join(base_dir, "RealignParameter", sub_id)
            )
        else:
            print(f"Skipping {sub_id}: Missing functional NIFTI or rp_*.txt file.")