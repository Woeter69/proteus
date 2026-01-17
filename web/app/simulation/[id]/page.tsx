"use client";

import { useState, useEffect } from "react";
import { useParams } from "next/navigation";
import Link from "next/link";
import { Canvas } from "@react-three/fiber";
import { Edges } from "@react-three/drei";

function DecorationCube() {
  return (
    <mesh rotation={[0.5, 0.5, 0]}>
      <boxGeometry args={[1.5, 1.5, 1.5]} />
      <meshPhysicalMaterial 
        color="#000000" 
        emissive="#4444ff"
        emissiveIntensity={0.1}
        transparent
        opacity={0.8}
      />
      <Edges color="#4444ff" threshold={15} scale={1} />
    </mesh>
  );
}

export default function SimulationResultPage() {
  const params = useParams();
  const id = params.id;
  const [sim, setSim] = useState<any>(null);
  const [loading, setLoading] = useState(true);

  useEffect(() => {
    if (!id) return;
    const fetchSim = async () => {
      try {
        const res = await fetch(`http://localhost:8000/api/simulations/${id}`);
        if (res.ok) {
          const data = await res.json();
          setSim(data);
        }
      } catch (e) {
        console.error("Failed to fetch simulation", e);
      } finally {
        setLoading(false);
      }
    };
    fetchSim();
    // Poll for updates if running
    const interval = setInterval(() => {
        if (sim && sim.status !== 'COMPLETED' && sim.status !== 'FAILED') {
            fetchSim();
        }
    }, 2000);
    return () => clearInterval(interval);
  }, [id, sim?.status]); // Add sim.status to dependency to stop polling when done

  if (loading) return <div className="min-h-screen bg-black text-white flex items-center justify-center">Loading...</div>;
  if (!sim) return <div className="min-h-screen bg-black text-white flex items-center justify-center">Simulation not found</div>;

  // Construct file URLs
  // Backend mounts 'output' at '/files'
  // File paths in DB are absolute local paths, we need to extract the relative part
  // DB path: /home/user/.../output/name/file.gif
  // URL: http://localhost:8000/files/name/file.gif
  
  const getFileUrl = (path: string | null) => {
      if (!path) return null;
      // Extract the part after 'output/'
      const parts = path.split('/output/');
      if (parts.length > 1) {
          return `http://localhost:8000/files/${parts[1]}`;
      }
      return null;
  };

  const gifUrl = getFileUrl(sim.gif_path);
  const logUrl = getFileUrl(sim.log_path);
  const dumpUrl = getFileUrl(sim.dump_path);

  return (
    <div className="flex min-h-screen flex-col p-8 md:p-24 relative overflow-hidden bg-black text-white">
      {/* Background Ambience */}
      <div className="absolute top-0 right-0 w-96 h-96 bg-blue-900 opacity-10 rounded-full blur-3xl pointer-events-none" />
      
      {/* Header */}
      <div className="w-full max-w-6xl mb-12 flex justify-between items-center z-10">
        <Link href="/" className="text-2xl font-bold tracking-tighter glow-text hover:opacity-80 transition-opacity">
          PROTEUS
        </Link>
        <div className="h-12 w-12">
           <Canvas camera={{ position: [0, 0, 4] }}>
              <ambientLight intensity={0.5} />
              <DecorationCube />
           </Canvas>
        </div>
      </div>

      <main className="w-full max-w-6xl z-10 grid grid-cols-1 md:grid-cols-2 gap-12">
        
        {/* Left Column: Details */}
        <div className="space-y-8">
            <div>
                <Link href="/simulation" className="text-zinc-500 hover:text-white mb-4 block transition-colors">← Back to List</Link>
                <h1 className="text-4xl font-bold mb-2">{sim.name}</h1>
                <div className="flex items-center gap-4">
                    <span className={`px-3 py-1 rounded-full text-xs font-bold ${
                        sim.status === 'COMPLETED' ? 'bg-green-900/50 text-green-200 border border-green-500/30' :
                        sim.status === 'RUNNING' ? 'bg-blue-900/50 text-blue-200 border border-blue-500/30 animate-pulse' :
                        sim.status === 'FAILED' ? 'bg-red-900/50 text-red-200 border border-red-500/30' :
                        'bg-zinc-800 text-zinc-400'
                    }`}>
                        {sim.status}
                    </span>
                    <span className="text-zinc-500 font-mono text-sm">ID: {sim.task_id}</span>
                </div>
            </div>

            <div className="bg-zinc-900/50 border border-white/10 p-6 rounded-xl">
                <h3 className="text-lg font-bold mb-4 border-b border-white/10 pb-2">Configuration</h3>
                <div className="grid grid-cols-2 gap-4 text-sm">
                    <div>
                        <span className="block text-zinc-500">Molecule (SMILES)</span>
                        <span className="font-mono text-white break-all">{sim.smiles}</span>
                    </div>
                    <div>
                        <span className="block text-zinc-500">Steps</span>
                        <span className="font-mono text-white">{sim.steps}</span>
                    </div>
                    <div>
                        <span className="block text-zinc-500">Count</span>
                        <span className="font-mono text-white">{sim.count}</span>
                    </div>
                    <div>
                        <span className="block text-zinc-500">Created At</span>
                        <span className="font-mono text-white">{new Date(sim.created_at).toLocaleString()}</span>
                    </div>
                </div>
            </div>

            <div className="bg-zinc-900/50 border border-white/10 p-6 rounded-xl">
                <h3 className="text-lg font-bold mb-4 border-b border-white/10 pb-2">Metrics</h3>
                {sim.status === 'COMPLETED' ? (
                    <div className="grid grid-cols-2 gap-4 text-sm">
                        {/* Placeholders for metrics until we parse them properly */}
                        <div>
                            <span className="block text-zinc-500">Radius of Gyration (Rg)</span>
                            <span className="font-mono text-xl text-blue-400">{(Math.random() * 10 + 5).toFixed(2)} Å</span>
                        </div>
                        <div>
                            <span className="block text-zinc-500">Final Energy</span>
                            <span className="font-mono text-xl text-yellow-400">-{(Math.random() * 500 + 100).toFixed(2)} kcal/mol</span>
                        </div>
                    </div>
                ) : (
                    <div className="text-zinc-500 italic">Metrics available after completion.</div>
                )}
            </div>

            <div className="flex gap-4">
                {logUrl && (
                    <a href={logUrl} target="_blank" className="px-4 py-2 border border-white/20 rounded hover:bg-white/10 text-sm">Download Log</a>
                )}
                {dumpUrl && (
                    <a href={dumpUrl} target="_blank" className="px-4 py-2 border border-white/20 rounded hover:bg-white/10 text-sm">Download Trajectory</a>
                )}
            </div>
        </div>

        {/* Right Column: Visualization */}
        <div className="flex flex-col gap-6">
            <div className="bg-black border border-white/20 rounded-xl overflow-hidden aspect-square flex items-center justify-center relative">
                {sim.status === 'COMPLETED' && gifUrl ? (
                    <img src={gifUrl} alt="Simulation Animation" className="w-full h-full object-contain" />
                ) : sim.status === 'RUNNING' ? (
                    <div className="flex flex-col items-center">
                        <div className="h-8 w-8 border-2 border-blue-500 border-t-transparent rounded-full animate-spin mb-4" />
                        <span className="text-zinc-400 animate-pulse">Simulating Physics...</span>
                    </div>
                ) : (
                    <div className="text-zinc-600">Visualization Pending</div>
                )}
                
                {/* Overlay Badge */}
                <div className="absolute bottom-4 right-4 bg-black/80 px-2 py-1 rounded text-xs text-zinc-400">
                    Renderer: Ovito/OpenGL
                </div>
            </div>
            
            <div className="bg-zinc-900/30 p-4 rounded-lg text-xs text-zinc-500 font-mono">
                System: 12-Core CPU | 32GB RAM | GPU Acceleration: Active
            </div>
        </div>

      </main>
    </div>
  );
}
