"use client";

import { Canvas, useFrame } from "@react-three/fiber";
import { Edges, Float, Stars, Sparkles, MeshTransmissionMaterial, Environment, Lightformer } from "@react-three/drei";
import { Bloom, EffectComposer } from "@react-three/postprocessing";
import { useRef } from "react";
import * as THREE from "three";
import Link from "next/link";

function FancyCube() {
  const groupRef = useRef<THREE.Group>(null);
  const coreRef = useRef<THREE.Group>(null);
  const ringRef = useRef<THREE.Group>(null);

  useFrame((state, delta) => {
    const t = state.clock.getElapsedTime();
    if (groupRef.current) {
      groupRef.current.rotation.y = Math.sin(t * 0.1) * 0.2;
      groupRef.current.rotation.z = Math.cos(t * 0.1) * 0.1;
    }
    if (coreRef.current) {
      coreRef.current.rotation.x = t * 0.5;
      coreRef.current.rotation.y = t * 0.7;
    }
    if (ringRef.current) {
      ringRef.current.rotation.x = t * 0.2;
      ringRef.current.rotation.y = t * 0.2;
    }
  });

  return (
    <group position={[2.5, 0, 0]} ref={groupRef}>
      <Float speed={4} rotationIntensity={0.5} floatIntensity={1}>
        
        {/* Outer Glass Cube */}
        <mesh>
          <boxGeometry args={[2, 2, 2]} />
          <MeshTransmissionMaterial 
            backside
            samples={16}
            resolution={1024}
            transmission={1}
            roughness={0.05}
            thickness={2.5}
            ior={1.7}
            chromaticAberration={2.0}
            anisotropy={0.5}
            distortion={0.5}
            distortionScale={0.5}
            temporalDistortion={0.1}
            clearcoat={1}
            attenuationDistance={5}
            attenuationColor="#ffffff"
            color="#eefbff"
          />
          <Edges threshold={1} scale={1.001}>
            <meshBasicMaterial color="white" transparent opacity={0.2} />
          </Edges>
        </mesh>

        {/* Inner Reactor Core */}
        <group ref={coreRef}>
            {/* Glowing Nucleus */}
            <mesh>
                <icosahedronGeometry args={[0.5, 1]} />
                <meshBasicMaterial color={[0, 10, 20]} toneMapped={false} wireframe />
            </mesh>
            <mesh>
                <icosahedronGeometry args={[0.4, 0]} />
                <meshBasicMaterial color={[0.5, 2, 5]} transparent opacity={0.2} />
            </mesh>
        </group>

        {/* Orbiting Rings */}
        <group ref={ringRef}>
            <mesh rotation={[Math.PI / 2, 0, 0]}>
                <torusGeometry args={[0.8, 0.015, 16, 64]} />
                <meshBasicMaterial color={[2, 10, 20]} toneMapped={false} />
            </mesh>
             <mesh rotation={[0, Math.PI / 2, 0]}>
                <torusGeometry args={[0.7, 0.015, 16, 64]} />
                <meshBasicMaterial color={[2, 20, 10]} toneMapped={false} />
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
          <ambientLight intensity={0.2} />
          <pointLight position={[20, 20, 20]} intensity={0.5} color="#ffffff" />
          
          <Stars radius={300} depth={60} count={20000} factor={8} saturation={0} fade speed={1} />
          
          <Environment resolution={1024}>
            <color attach="background" args={["#000000"]} />
            {/* BAKING STARS INTO ENVIRONMENT FOR REFRACTION */}
            <group rotation={[0, 0, 0]}>
               <Stars radius={100} depth={50} count={5000} factor={4} saturation={0} fade speed={1} />
               {/* Bright Sparkles to create sharp Chromatic Aberration points */}
               <Sparkles count={1000} scale={150} size={15} speed={0} opacity={1} color="#ffffff" />
               <Sparkles count={500} scale={100} size={25} speed={0} opacity={1} color="#ffeebb" />
            </group>
          </Environment>
          
          <FancyCube />

          <EffectComposer enableNormalPass={false}>
            <Bloom 
              luminanceThreshold={0.5} 
              mipmapBlur 
              intensity={2.0} 
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