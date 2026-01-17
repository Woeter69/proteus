"use client";

import { useState, useEffect } from "react";
import { Canvas } from "@react-three/fiber";
import { Edges } from "@react-three/drei";
import Link from "next/link";

function DecorationCube() {
  return (
    <mesh rotation={[0.5, 0.5, 0]}>
      <boxGeometry args={[1.5, 1.5, 1.5]} />
      <meshPhysicalMaterial 
        color="#000000" 
        emissive="#ffffff"
        emissiveIntensity={0.2}
        transparent
        opacity={0.8}
      />
      <Edges color="#ffffff" threshold={15} scale={1} />
    </mesh>
  );
}

function RecentSimulations() {
  const [sims, setSims] = useState<any[]>([]);

  const fetchSims = async () => {
    try {
      const res = await fetch("http://localhost:8000/api/simulations");
      if (res.ok) {
        const data = await res.json();
        setSims(data);
      }
    } catch (e) {
      console.error("Failed to fetch simulations", e);
    }
  };

  useEffect(() => {
    fetchSims();
    const interval = setInterval(fetchSims, 2000);
    return () => clearInterval(interval);
  }, []);

  return (
    <div className="w-full max-w-4xl mt-12">
      <h2 className="text-2xl font-bold mb-4 text-white">Recent Simulations</h2>
      <div className="grid gap-4">
        {sims.map((sim) => (
          <div key={sim.id} className="bg-zinc-900/50 border border-white/10 p-4 rounded-lg flex items-center justify-between hover:bg-zinc-800/50 transition-all">
            <div>
              <div className="font-bold text-lg text-white">{sim.name}</div>
              <div className="text-sm text-zinc-400 font-mono">{sim.smiles}</div>
              <div className="text-xs text-zinc-500 mt-1">ID: {sim.task_id}</div>
            </div>
            <div className="flex flex-col items-end gap-2">
              <span className={`px-3 py-1 rounded-full text-xs font-bold ${
                sim.status === 'COMPLETED' ? 'bg-green-900/50 text-green-200 border border-green-500/30' :
                sim.status === 'RUNNING' ? 'bg-blue-900/50 text-blue-200 border border-blue-500/30 animate-pulse' :
                sim.status === 'FAILED' ? 'bg-red-900/50 text-red-200 border border-red-500/30' :
                'bg-zinc-800 text-zinc-400'
              }`}>
                {sim.status}
              </span>
              {sim.status === 'COMPLETED' && (
                <Link href={`/simulation/${sim.id}`} className="text-xs text-blue-400 hover:text-blue-300 underline">
                  View Results
                </Link>
              )}
            </div>
          </div>
        ))}
        {sims.length === 0 && (
          <div className="text-zinc-500 text-center py-8">No simulations found. Start one above!</div>
        )}
      </div>
    </div>
  );
}

export default function SimulationPage() {
  const [smiles, setSmiles] = useState("");
  const [name, setName] = useState("");
  const [isSubmitting, setIsSubmitting] = useState(false);
  const [taskId, setTaskId] = useState<string | null>(null);
  const [error, setError] = useState<string | null>(null);

  const handleSubmit = async (e: React.FormEvent) => {
    e.preventDefault();
    setIsSubmitting(true);
    setError(null);
    setTaskId(null);
    
    try {
      const response = await fetch("http://localhost:8000/api/simulate", {
        method: "POST",
        headers: {
          "Content-Type": "application/json",
        },
        body: JSON.stringify({
          smiles,
          name: name || `sim_${Date.now()}`, 
          render: true 
        }),
      });

      if (!response.ok) {
        throw new Error("Failed to start simulation");
      }

      const data = await response.json();
      setTaskId(data.task_id);
      // Clear form on success
      setSmiles("");
      setName("");
    } catch (err) {
      setError(err instanceof Error ? err.message : "An unknown error occurred");
    } finally {
      setIsSubmitting(false);
    }
  };

  return (
    <div className="flex min-h-screen flex-col items-center p-8 md:p-24 relative overflow-hidden bg-black text-white">
      {/* Background Ambience */}
      <div className="absolute top-0 right-0 w-96 h-96 bg-white opacity-5 rounded-full blur-3xl pointer-events-none" />
      <div className="absolute bottom-0 left-0 w-64 h-64 bg-white opacity-5 rounded-full blur-3xl pointer-events-none" />

      {/* Header */}
      <div className="w-full max-w-4xl mb-12 flex justify-between items-center z-10">
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

      {/* Form Container */}
      <main className="w-full max-w-xl z-10 mb-12">
        <div className="bg-zinc-900/50 border border-white/10 p-8 rounded-2xl backdrop-blur-md shadow-2xl">
          <h1 className="text-3xl font-bold mb-2">New Simulation</h1>
          <p className="text-zinc-400 mb-8">Enter a generic SMILES string to generate a polymer topology.</p>

          <form onSubmit={handleSubmit} className="space-y-6">
            <div className="space-y-2">
              <label htmlFor="name" className="block text-sm font-medium text-zinc-300">
                Simulation Name (Optional)
              </label>
              <input
                id="name"
                type="text"
                value={name}
                onChange={(e) => setName(e.target.value)}
                placeholder="e.g., Polyethylene Test"
                className="w-full bg-black/50 border border-white/20 rounded-lg px-4 py-3 text-white placeholder-zinc-600 focus:outline-none focus:ring-2 focus:ring-white/50 focus:border-white/50 transition-all"
              />
            </div>

            <div className="space-y-2">
              <label htmlFor="smiles" className="block text-sm font-medium text-zinc-300">
                SMILES String
              </label>
              <input
                id="smiles"
                type="text"
                value={smiles}
                onChange={(e) => setSmiles(e.target.value)}
                placeholder="e.g., C=CC"
                className="w-full bg-black/50 border border-white/20 rounded-lg px-4 py-3 text-white placeholder-zinc-600 focus:outline-none focus:ring-2 focus:ring-white/50 focus:border-white/50 transition-all font-mono"
                required
              />
            </div>

            <div className="pt-4">
              <button
                type="submit"
                disabled={isSubmitting}
                className="w-full py-4 bg-white text-black font-bold rounded-lg hover:bg-zinc-200 transition-all disabled:opacity-50 disabled:cursor-not-allowed flex items-center justify-center gap-2 glow-box"
              >
                {isSubmitting ? (
                  <>
                    <div className="h-4 w-4 border-2 border-black border-t-transparent rounded-full animate-spin" />
                    Processing...
                  </>
                ) : (
                  "Initialize Simulation"
                )}
              </button>
            </div>
            
            {error && (
              <div className="p-4 bg-red-900/20 border border-red-500/50 rounded-lg text-red-200 text-sm">
                Error: {error}
              </div>
            )}

            {taskId && (
              <div className="p-4 bg-green-900/20 border border-green-500/50 rounded-lg text-green-200 text-sm animate-pulse">
                Simulation Started! Check the list below for updates.
              </div>
            )}

            <div className="text-center">
               <p className="text-xs text-zinc-500 mt-4">
                 Standard simulation runs on CPU. GPU acceleration enabled if available.
               </p>
            </div>
          </form>
        </div>
      </main>

      {/* Recent Simulations List */}
      <RecentSimulations />
    </div>
  );
}
