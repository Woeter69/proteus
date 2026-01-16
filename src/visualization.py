"""
Module V: Visualization
Uses Ovito to render the simulation trajectory into an animation.
"""

import sys
from pathlib import Path
import warnings

# Suppress OVITO PyPI warning as requested by the library
warnings.filterwarnings('ignore', message='.*OVITO.*PyPI')

# Debugging Import
try:
    from ovito.io import import_file
    from ovito.vis import Viewport, TachyonRenderer, OpenGLRenderer
    from ovito.modifiers import AmbientOcclusionModifier, ColorCodingModifier
    OVITO_AVAILABLE = True
    IMPORT_ERROR = None
except ImportError as e:
    OVITO_AVAILABLE = False
    IMPORT_ERROR = e

def render_trajectory(dump_path: Path, output_gif: Path):
    """
    Renders the trajectory.dump file into a GIF animation using Ovito.
    """
    print(f"[*] Visualization Module Invoked for: {dump_path}")
    
    if not OVITO_AVAILABLE:
        print(f"Warning: Ovito python module not found. Skipping visualization.")
        return

    if not dump_path.exists():
        print(f"Error: Dump file not found at {dump_path}")
        return

    print(f"[*] Rendering visualization to {output_gif.name}...")

    try:
        # 1. Load Data
        pipeline = import_file(str(dump_path))
        pipeline.add_to_scene()
        
        num_frames = pipeline.source.num_frames
        print(f"[*] Loaded Trajectory: {num_frames} frames found.")

        if num_frames == 0:
            print("Warning: Trajectory has 0 frames. Nothing to render.")
            return

        # 2. visual settings
        # Ambient Occlusion is expensive on CPU, but fast on GPU if supported
        pipeline.modifiers.append(AmbientOcclusionModifier(intensity=0.3))
        
        try:
            pipeline.modifiers.append(ColorCodingModifier(property='Molecule Identifier'))
        except Exception:
            pass

        # 3. Setup Camera/Viewport
        vp = Viewport()
        vp.type = Viewport.Type.Perspective
        vp.zoom_all()
        
        # 4. Render Animation
        print("[*] Starting render... (Trying GPU/OpenGL acceleration)")
        
        try:
            # Try GPU Renderer first
            vp.render_anim(
                filename=str(output_gif),
                size=(800, 600),
                fps=10,
                renderer=OpenGLRenderer()
            )
            print("[*] GPU Rendering finished.")
            
        except Exception as e:
            print(f"Warning: GPU Rendering failed ({e}). Falling back to CPU Raytracing...")
            # Fallback to CPU Tachyon
            vp.render_anim(
                filename=str(output_gif),
                size=(800, 600),
                fps=10,
                renderer=TachyonRenderer(antialiasing_samples=1)
            )
            print("[*] CPU Rendering finished.")
        
        print(f"[*] Visualization saved successfully: {output_gif}")
        print("[*] Rendering finished.")
        
        print(f"[*] Visualization saved successfully: {output_gif}")
        
    except Exception as e:
        print(f"Error rendering visualization: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    pass