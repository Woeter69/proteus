"use client";

import { Canvas, useFrame } from "@react-three/fiber";
import { Float, Stars, Environment, Edges, MeshTransmissionMaterial } from "@react-three/drei";
import { Bloom, EffectComposer } from "@react-three/postprocessing";
import { useRef } from "react";
import * as THREE from "three";
import Link from "next/link";

function IridescentReactor() {
  const outerRef = useRef<THREE.Mesh>(null);

  useFrame((state) => {
    const t = state.clock.getElapsedTime();
    if (outerRef.current) {
      outerRef.current.rotation.x = t * 0.1;
      outerRef.current.rotation.y = t * 0.15;
    }
  });

  return (
    <group position={[2.5, 0, 0]}>
      <Float speed={2} rotationIntensity={0.2} floatIntensity={0.5}>
        
        {/* Glass Cube with Chromatic Aberration */}
        <mesh ref={outerRef}>
          <boxGeometry args={[2.2, 2.2, 2.2]} />
          <MeshTransmissionMaterial
            backside
            samples={16}
            resolution={1024}
            transmission={1}
            thickness={1.5}
            roughness={0.05}
            chromaticAberration={0.8}
            anisotropy={0.1}
            distortion={1.2}
            distortionScale={1.0}
            temporalDistortion={0.4}
            ior={1.7}
            color="#ffffff"
            attenuationDistance={1}
            attenuationColor="#ffffff"
          />
          <Edges color="white" scale={1.02} threshold={1} />
        </mesh>

      </Float>
    </group>
  );
}

function ManualStars() {
  // Generate random positions for stars
  const stars = new Array(200).fill(0).map(() => ({
    pos: [
      (Math.random() - 0.5) * 20,
      (Math.random() - 0.5) * 20,
      (Math.random() - 0.5) * 20
    ] as [number, number, number],
    scale: Math.random() * 0.05 + 0.02
  }));

  return (
    <group>
      {stars.map((s, i) => (
        <mesh key={i} position={s.pos}>
          <sphereGeometry args={[s.scale, 8, 8]} />
          <meshBasicMaterial color="white" toneMapped={false} />
        </mesh>
      ))}
    </group>
  );
}

export default function Home() {
  return (
    <div className="flex min-h-screen w-full flex-col md:flex-row items-center justify-center p-8 md:p-24 relative overflow-hidden bg-black">
      
      {/* Full-screen 3D Background */}
      <div className="absolute inset-0 z-0">
        <Canvas 
          shadows
          dpr={[1, 2]}
          camera={{ position: [0, 0, 8], fov: 40 }}
        >
          <ambientLight intensity={0.5} />
          
          {/* Background Stars (Visual) */}
          <Stars radius={100} depth={50} count={5000} factor={6} saturation={0} fade speed={1} />

          {/* Reflection Environment (Actual bright objects) */}
          <Environment resolution={1024}>
            <ManualStars />
            <color attach="background" args={["#000000"]} />
          </Environment>
          
          <IridescentReactor />

          <EffectComposer enableNormalPass={false}>
            <Bloom 
              luminanceThreshold={0.8} 
              mipmapBlur 
              intensity={0.5} 
              radius={0.4} 
            />
          </EffectComposer>
        </Canvas>
      </div>

      {/* Foreground Content */}
      <div className="flex-1 z-10 flex flex-col items-center md:items-start text-center md:text-left space-y-6">
        <h1 className="text-6xl md:text-8xl font-bold tracking-tighter glow-text text-white">
          PROTEUS
        </h1>
        <p className="text-xl md:text-2xl text-gray-300 max-w-md font-light">
          Automated Molecular Dynamics Pipeline for Polymer Engineering.
        </p>
        <div className="flex gap-4 pt-4">
          <Link href="/simulation" className="px-8 py-3 bg-white text-black font-bold rounded-full hover:bg-gray-200 transition-all glow-box">
            Start Simulation
          </Link>
          <button className="px-8 py-3 border border-white text-white font-bold rounded-full hover:bg-white/10 transition-all">
            Documentation
          </button>
        </div>
      </div>

      {/* Spacer for the Cube which is now positioned via Three.js group on the right */}
      <div className="flex-1 pointer-events-none" />
      
    </div>
  );
}
