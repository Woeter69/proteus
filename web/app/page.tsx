"use client";

import { Canvas, useFrame } from "@react-three/fiber";
import { Edges, Float, Stars, Sparkles, MeshTransmissionMaterial, Environment } from "@react-three/drei";
import { Bloom, EffectComposer } from "@react-three/postprocessing";
import { useRef } from "react";
import * as THREE from "three";
import Link from "next/link";

function FancyCube() {
  const spinnerRef = useRef<THREE.Group>(null);
  const innerRef = useRef<THREE.Group>(null);
  const coreRef = useRef<THREE.Group>(null);

  useFrame((state, delta) => {
    if (spinnerRef.current) {
      spinnerRef.current.rotation.y -= delta * 0.4;
    }
    if (innerRef.current) {
      innerRef.current.rotation.y += delta * 0.8;
      innerRef.current.rotation.z += delta * 0.4;
    }
    if (coreRef.current) {
      coreRef.current.rotation.x -= delta * 1.2;
      coreRef.current.rotation.y -= delta * 0.6;
    }
  });

  const cubeRotation: [number, number, number] = [Math.atan(1 / Math.sqrt(2)), 0, Math.PI / 4];
  const size = 1.8;
  const innerSize = size * 0.5;
  const coreSize = size * 0.2;

  return (
    <group position={[2.5, 0, 0]}> 
      <Float speed={2} rotationIntensity={0.2} floatIntensity={0.5}>
        <group ref={spinnerRef}>
           
           {/* Outer Cube (Prism Glass) */}
           <mesh rotation={cubeRotation}>
              <boxGeometry args={[size, size, size]} />
              <MeshTransmissionMaterial
                backside
                samples={32}
                resolution={2048}
                transmission={1}
                roughness={0.0}
                thickness={0.2} 
                ior={1.5}
                chromaticAberration={1.0}
                anisotropy={0.1} 
                distortion={0.0} // Removed distortion to see inside
                distortionScale={0.0}
                temporalDistortion={0.0} 
                clearcoat={1}
                attenuationDistance={5}
                attenuationColor="#ffffff"
                color="#ffffff"
              />
              <Edges threshold={1} scale={1.001}>
                 <meshBasicMaterial color={[3, 3, 3]} toneMapped={false} />
              </Edges>
           </mesh>

           {/* Middle Cube (Glowing Outlines Only) */}
           <group ref={innerRef} rotation={cubeRotation}>
             <mesh>
               <boxGeometry args={[innerSize, innerSize, innerSize]} />
               <meshBasicMaterial color="#000000" transparent opacity={0} /> {/* Invisible faces */}
               <Edges threshold={1} scale={1.001}>
                 <meshBasicMaterial color={[10, 10, 10]} toneMapped={false} />
               </Edges>
             </mesh>
           </group>

           {/* Core Cube (Small Glowing Outline) */}
           <group ref={coreRef} rotation={cubeRotation}>
             <mesh>
               <boxGeometry args={[coreSize * 1.5, coreSize * 1.5, coreSize * 1.5]} />
               <meshBasicMaterial color="#000000" transparent opacity={0} />
               <Edges threshold={1} scale={1.001}>
                 <meshBasicMaterial color={[15, 15, 15]} toneMapped={false} />
               </Edges>
             </mesh>
           </group>

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
          <pointLight position={[10, 10, 10]} intensity={1} color="#ffffff" />
          
          <Stars radius={300} depth={60} count={20000} factor={8} saturation={0} fade speed={1} />
          
          <Environment resolution={1024}>
            <group rotation={[-Math.PI / 4, 0, 0]}>
              <mesh scale={200}>
                <sphereGeometry args={[1, 64, 64]} />
                <meshBasicMaterial color="#000000" side={THREE.BackSide} />
              </mesh>
              <Sparkles count={500} scale={50} size={6} speed={0} opacity={1} color="#ffffff" />
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