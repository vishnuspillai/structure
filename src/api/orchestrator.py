import subprocess
import asyncio
import sys
import os

class PipelineOrchestrator:
    def __init__(self, root_dir):
        self.root_dir = root_dir
        self.steps = [
            ("Module 1: Extract and filter rare baseline variants", "src/module1_variant_mining/variant_mining.py"),
            ("Module 2: Phase A/B/C Annotations", "src/module2_annotation/annotation_layer.py"),
            ("Module 3: Structural Mapping", "src/module3_spatial/spatial_annotation.py"),
            ("Module 4: Scorer/Prioritization", "src/module4_prioritization/prioritization.py"),
            ("Module 5: Mechanistic Deep Dive", "src/module5_mechanistic/mechanistic_metrics.py"),
            ("Module 6: Clinical Review", "src/module6_clinical/clinical_review.py"),
            ("Module 7: CI Validation", "src/module7_validation/structural_enrichment_ci.py")
        ]
        self.current_process = None
        self.output_queue = asyncio.Queue()

    async def run_step(self, step_index):
        if step_index < 0 or step_index >= len(self.steps):
            return False, "Invalid step index"
        
        desc, script_path = self.steps[step_index]
        abs_script_path = os.path.join(self.root_dir, script_path)
        
        process = await asyncio.create_subprocess_exec(
            sys.executable, abs_script_path,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            cwd=self.root_dir
        )
        self.current_process = process
        
        async for line in process.stdout:
            await self.output_queue.put(line.decode('utf-8', errors='replace').strip())
            
        await process.wait()
        return process.returncode == 0, desc

    def get_steps(self):
        return [{"index": i, "name": name} for i, (name, _) in enumerate(self.steps)]
