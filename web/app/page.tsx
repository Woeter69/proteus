"use client";

import { Canvas, useFrame } from "@react-three/fiber";
import { Edges, Float, Stars, Sparkles, MeshTransmissionMaterial, Environment, Line } from "@react-three/drei";
import { Bloom, EffectComposer } from "@react-three/postprocessing";
import { useRef, useMemo } from "react";
import * as THREE from "three";
import Link from "next/link";

function FancyCube() {
  const spinnerRef = useRef<THREE.Group>(null);

  useFrame((state, delta) => {
    if (spinnerRef.current) {
      spinnerRef.current.rotation.y -= delta * 0.4;
    }
  });

  // Calculate rotation to align the cube's diagonal with the Y-axis (isometric look)
  const cubeRotation: [number, number, number] = [Math.atan(1 / Math.sqrt(2)), 0, Math.PI / 4];
  const size = 1.8;
  const halfSize = size / 2;
  const trapOffset = 0.1; // Reduced from 0.12

  // Define the trapezium shape that expands from the cube edges onto the faces
  const trapShape = useMemo(() => {
    const s = new THREE.Shape();
    s.moveTo(-halfSize, halfSize);
    s.lineTo(halfSize, halfSize);
    s.lineTo(halfSize * 0.7, halfSize - trapOffset);
    s.lineTo(-halfSize * 0.7, halfSize - trapOffset);
    s.closePath();
    return s;
  }, [halfSize, trapOffset]);

  // Define the twin blade sword shape (a long, thin diamond)
  const bladeShape = useMemo(() => {
    const s = new THREE.Shape();
    s.moveTo(0, halfSize);
    s.lineTo(0.03, 0); // Narrower width
    s.lineTo(0, -halfSize);
    s.lineTo(-0.03, 0); // Narrower width
    s.closePath();
    return s;
  }, [halfSize]);

  return (
    <group position={[2.5, 0, 0]}> {/* Position the cube to the right */}
      <Float speed={2} rotationIntensity={0.2} floatIntensity={0.5}>
        <group ref={spinnerRef}>
           <mesh rotation={cubeRotation}>
              <boxGeometry args={[size, size, size]} />
              <MeshTransmissionMaterial
                backside
                samples={32}
                resolution={2048}
                transmission={1}
                roughness={0.0}
                thickness={1.5}
                ior={2.4}
                chromaticAberration={0.6}
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
           
           <group rotation={cubeRotation}>
              <Edges color="#ffffff" threshold={15} scale={1} />

              {[0, 1, 2, 3, 4, 5].map((i) => {
                const rot: [number, number, number] = [0, 0, 0];
                if (i === 4) rot[0] = Math.PI / 2;
                else if (i === 5) rot[0] = -Math.PI / 2;
                else rot[1] = (i * Math.PI) / 2;

                return (
                  <group key={`face-${i}`} rotation={rot}>
                    <group position={[0, 0, halfSize * 1.01]}>
                      <mesh>
                        <planeGeometry args={[0.4, 0.4]} />
                        <meshBasicMaterial color="#ffffff" transparent opacity={0.9} />
                      </mesh>

                      {[0, 1].map((r) => (
                        <group key={`face-blade-${r}`} rotation={[0, 0, r * Math.PI / 2]}>
                           <mesh>
                             <shapeGeometry args={[bladeShape]} />
                             <meshBasicMaterial color="#ffffff" transparent opacity={0.6} />
                           </mesh>
                           <Line
                             points={[[0, -halfSize, 0], [0, halfSize, 0]]}
                             color="#ffffff"
                             lineWidth={1.5}
                             transparent
                             opacity={0.8}
                           />
                        </group>
                      ))}
                      
                      {[0, 1, 2, 3].map((side) => (
                        <group key={`trap-${side}`} rotation={[0, 0, (side * Math.PI) / 2]}>
                           <mesh>
                             <shapeGeometry args={[trapShape]} />
                             <meshBasicMaterial color="#ffffff" transparent opacity={0.4} />
                           </mesh>
                           <Line
                             points={[
                               [-halfSize, halfSize, 0],
                               [halfSize, halfSize, 0],
                               [halfSize * 0.7, halfSize - trapOffset, 0],
                               [-halfSize * 0.7, halfSize - trapOffset, 0],
                               [-halfSize, halfSize, 0],
                             ]}
                             color="#ffffff"
                             lineWidth={1}
                             transparent
                             opacity={0.8}
                           />
                        </group>
                      ))}
                    </group>
                  </group>
                );
              })}
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
          <Sparkles scale={100} count={1000} size={2} speed={0.2} opacity={0.15} color="#ffffff" />
          <Sparkles scale={50} count={500} size={1} speed={0.1} opacity={0.3} color="#ffffff" />
          <Sparkles scale={20} count={200} size={1.5} speed={0.3} opacity={0.5} color="#ffffff" />
          
          <Environment resolution={1024}>
            <group rotation={[-Math.PI / 4, 0, 0]}>
              <mesh scale={200}>
                <sphereGeometry args={[1, 64, 64]} />
                <meshBasicMaterial color="#000000" side={THREE.BackSide} />
              </mesh>
              <Sparkles count={2000} scale={100} size={4} speed={0} opacity={1} color="#ffffff" />
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
