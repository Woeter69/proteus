"use client";

import { Canvas, useFrame } from "@react-three/fiber";
import { Float, Stars, Environment, Lightformer } from "@react-three/drei";
import { Bloom, EffectComposer } from "@react-three/postprocessing";
import { useRef } from "react";
import * as THREE from "three";
import Link from "next/link";

function IridescentReactor() {
  const outerRef = useRef<THREE.Mesh>(null);
  const coreRef = useRef<THREE.Group>(null);
  const ringRef = useRef<THREE.Group>(null);

  useFrame((state) => {
    const t = state.clock.getElapsedTime();
    if (outerRef.current) {
      outerRef.current.rotation.x = t * 0.1;
      outerRef.current.rotation.y = t * 0.15;
    }
    if (coreRef.current) {
      coreRef.current.rotation.x = -t * 0.2;
      coreRef.current.rotation.z = t * 0.1;
    }
    if (ringRef.current) {
      ringRef.current.rotation.y = t * 0.4;
      ringRef.current.rotation.x = Math.sin(t) * 0.2;
    }
  });

  return (
    <group position={[2.5, 0, 0]}>
      <Float speed={2} rotationIntensity={0.2} floatIntensity={0.5}>
        
        {/* Iridescent Outer Shell - Dark base for max rainbow contrast */}
        <mesh ref={outerRef}>
          <boxGeometry args={[2.2, 2.2, 2.2]} />
          <meshPhysicalMaterial
            transparent
            transmission={0.8}
            opacity={1}
            roughness={0.1}
            metalness={0.1}
            ior={1.5}
            iridescence={1}
            iridescenceIOR={1.8}
            iridescenceThicknessRange={[100, 800]}
            envMapIntensity={1.5}
            clearcoat={1}
            color="#101010" 
          />
        </mesh>

        {/* Inner Structure - Gold/Neutral to avoid "Blue Neon" look */}
        <group ref={coreRef}>
          <mesh>
            <icosahedronGeometry args={[0.7, 0]} />
            <meshStandardMaterial color="#ffccaa" wireframe emissive="#ff8800" emissiveIntensity={0.5} />
          </mesh>
          <mesh>
            <sphereGeometry args={[0.4, 32, 32]} />
            <meshBasicMaterial color="#ffffff" transparent opacity={0.2} />
          </mesh>
        </group>

        {/* Orbiting Rings - Metallic Silver */}
        <group ref={ringRef}>
          <mesh rotation={[Math.PI / 2, 0, 0]}>
            <torusGeometry args={[1.3, 0.02, 16, 100]} />
            <meshStandardMaterial color="#aaaaaa" metalness={1} roughness={0.2} />
          </mesh>
          <mesh rotation={[0, Math.PI / 2, 0]}>
            <torusGeometry args={[1.1, 0.02, 16, 100]} />
            <meshStandardMaterial color="#aaaaaa" metalness={1} roughness={0.2} />
          </mesh>
        </group>

      </Float>
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
          
          {/* Neutral white lighting */}
          <Environment resolution={1024}>
            <color attach="background" args={["#000000"]} />
            <group rotation={[0, 0, 0]}>
               <Stars radius={100} depth={50} count={5000} factor={4} saturation={0} fade speed={1} />
               
               {/* Broad White Lightformers to drive the Iridescence Gradient */}
               <Lightformer form="rect" intensity={5} position={[5, 5, 5]} scale={[10, 10, 1]} color="white" />
               <Lightformer form="rect" intensity={5} position={[-5, 5, 5]} scale={[10, 10, 1]} color="white" />
               <Lightformer form="rect" intensity={5} position={[0, -5, 5]} scale={[10, 10, 1]} color="white" />
               <Lightformer form="circle" intensity={2} position={[0, 0, -10]} scale={[20, 20, 1]} color="white" />
            </group>
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
