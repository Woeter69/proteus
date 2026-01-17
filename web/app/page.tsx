"use client";

import { Canvas, useFrame } from "@react-three/fiber";
import { Edges, Float, Stars, Sparkles, MeshTransmissionMaterial, Environment, Line, Lightformer } from "@react-three/drei";
import { Bloom, EffectComposer } from "@react-three/postprocessing";
import { useRef } from "react";
import * as THREE from "three";
import Link from "next/link";

function FancyCube() {
  const spinnerRef = useRef<THREE.Group>(null);

  useFrame((state, delta) => {
    if (spinnerRef.current) {
      spinnerRef.current.rotation.y -= delta * 0.4;
    }
  });

  const cubeRotation: [number, number, number] = [Math.atan(1 / Math.sqrt(2)), 0, Math.PI / 4];

  return (
    <group position={[0, 0, 0]}>
      <Float speed={2} rotationIntensity={0.2} floatIntensity={0.5}>
        <group ref={spinnerRef}>
           {/* The Glass Cube */}
           <mesh rotation={cubeRotation}>
              <boxGeometry args={[2.5, 2.5, 2.5]} />
              <MeshTransmissionMaterial
                backside
                samples={32}
                resolution={2048}
                transmission={1}
                roughness={0.0}
                thickness={1.5}
                ior={2.4} // Diamond-like brilliance
                chromaticAberration={0.6} // Realistic dispersion
                anisotropy={0.3}
                distortion={0}
                distortionScale={0}
                temporalDistortion={0}
                clearcoat={1}
                attenuationDistance={5}
                attenuationColor="#ffffff"
                color="#ffffff"
                //@ts-ignore
                iridescence={0.8}
                //@ts-ignore
                iridescenceIOR={1.3}
                //@ts-ignore
                iridescenceThicknessRange={[100, 800]}
              />
           </mesh>
           
           {/* Subtle Glowy Edge */}
           <mesh rotation={cubeRotation}>
              <boxGeometry args={[2.51, 2.51, 2.51]} />
              <meshBasicMaterial transparent opacity={0} />
              <Edges color="#ffffff" threshold={15} scale={1} />
           </mesh>
        </group>
      </Float>
      
      <Sparkles count={15} scale={8} size={1} speed={0.2} opacity={0.3} color="#ffffff" />
    </group>
  );
}

export default function Home() {
  return (
    <div className="flex min-h-screen w-full flex-col md:flex-row items-center justify-center p-8 md:p-24 relative overflow-hidden bg-black">
      
      {/* Background ambient glow bits */}
      <div className="absolute top-1/4 left-1/4 w-96 h-96 bg-white opacity-5 rounded-full blur-3xl pointer-events-none" />
      <div className="absolute bottom-1/4 right-1/4 w-64 h-64 bg-white opacity-5 rounded-full blur-3xl pointer-events-none" />

      {/* Left Side: Text */}
      <div className="flex-1 z-10 flex flex-col items-center md:items-start text-center md:text-left space-y-6">
        <h1 className="text-6xl md:text-8xl font-bold tracking-tighter glow-text text-white">
          PROTEUS
        </h1>
        <p className="text-xl md:text-2xl text-gray-300 max-w-md font-light">
          Advanced automated pipeline for polymer simulation and molecular dynamics.
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

      {/* Right Side: 3D Scene */}
      <div className="flex-1 h-[600px] w-full z-10">
        <Canvas 
          shadows
          dpr={[1, 2]}
          camera={{ position: [0, 2, 8], fov: 40 }}
        >
          <ambientLight intensity={0.2} />
          <pointLight position={[10, 10, 10]} intensity={1} color="#ffffff" />
          
          {/* Deep Background Stars */}
          <Stars radius={100} depth={100} count={8000} factor={4} saturation={0} fade speed={1} />
          {/* Foreground Sparkles for parallax depth */}
          <Sparkles scale={100} count={500} size={2} speed={0.2} opacity={0.5} color="#ffffff" />
          
          {/* Environment reflecting ONLY stars */}
          <Environment resolution={1024}>
            <group rotation={[-Math.PI / 4, 0, 0]}>
              {/* The "Void" */}
              <mesh scale={100}>
                <sphereGeometry args={[1, 64, 64]} />
                <meshBasicMaterial color="#000000" side={THREE.BackSide} />
              </mesh>
              
              {/* Dense Starfield for reflections ONLY */}
              <Sparkles count={1000} scale={50} size={4} speed={0} opacity={1} color="#ffffff" />
              <Sparkles count={500} scale={30} size={2} speed={0} opacity={0.8} color="#ffffff" />
            </group>
          </Environment>

          <FancyCube />

          <EffectComposer disableNormalPass>
            <Bloom 
              luminanceThreshold={1} 
              mipmapBlur 
              intensity={1.0} 
              radius={0.3} 
            />
          </EffectComposer>

        </Canvas>
      </div>
    </div>
  );
}
