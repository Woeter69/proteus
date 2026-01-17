"use client";

import { Canvas, useFrame } from "@react-three/fiber";
import { Stars, Sparkles, Environment, Float, Edges, Lightformer, MeshTransmissionMaterial } from "@react-three/drei";
import { Bloom, EffectComposer } from "@react-three/postprocessing";
import { useRef, useMemo, Suspense } from "react";
import * as THREE from "three";
import Link from "next/link";

// Custom Shader for the Additive Spectral Ghost Layer
const SpectralGhostShader = {
  uniforms: {
    uTime: { value: 0 },
  },
  vertexShader: `
    varying vec2 vUv;
    varying vec3 vNormal;
    varying vec3 vPosition;
    void main() {
      vUv = uv;
      vNormal = normalize(normalMatrix * normal);
      vPosition = position;
      gl_Position = projectionMatrix * modelViewMatrix * vec4(position, 1.0);
    }
  `,
  fragmentShader: `
    uniform float uTime;
    varying vec2 vUv;
    varying vec3 vNormal;
    varying vec3 vPosition;

    vec3 hsv2rgb(vec3 c) {
      vec4 K = vec4(1.0, 2.0 / 3.0, 1.0 / 3.0, 3.0);
      vec3 p = abs(fract(c.xxx + K.xyz) * 6.0 - K.www);
      return c.z * mix(K.xxx, clamp(p - K.xxx, 0.0, 1.0), c.y);
    }

    // Simplex 2D noise
    vec3 permute(vec3 x) { return mod(((x*34.0)+1.0)*x, 289.0); }
    float snoise(vec2 v) {
      const vec4 C = vec4(0.211324865405187, 0.366025403784439, -0.577350269189626, 0.024390243902439);
      vec2 i  = floor(v + dot(v, C.yy) );
      vec2 x0 = v -   i + dot(i, C.xx);
      vec2 i1 = (x0.x > x0.y) ? vec2(1.0, 0.0) : vec2(0.0, 1.0);
      vec4 x12 = x0.xyxy + C.xxzz;
      x12.xy -= i1;
      i = mod(i, 289.0);
      vec3 p = permute( permute( i.y + vec3(0.0, i1.y, 1.0 )) + i.x + vec3(0.0, i1.x, 1.0 ));
      vec3 m = max(0.5 - vec3(dot(x0,x0), dot(x12.xy,x12.xy), dot(x12.zw,x12.zw)), 0.0);
      m = m*m ; m = m*m ;
      vec3 x = 2.0 * fract(p * C.www) - 1.0;
      vec3 h = abs(x) - 0.5;
      vec3 ox = floor(x + 0.5);
      vec3 a0 = x - ox;
      m *= 1.79284291400159 - 0.85373472095314 * ( a0*a0 + h*h );
      vec3 g;
      g.x  = a0.x  * x0.x  + h.x  * x0.y;
      g.yz = a0.yz * x12.xz + h.yz * x12.yw;
      return 130.0 * dot(m, g);
    }

    void main() {
      vec2 uv = vUv;
      float noise = snoise(uv * 2.0 + uTime * 0.1);
      
      float hue = fract(uv.x + uv.y * 0.3 + noise * 0.2 + uTime * 0.05);
      vec3 rainbow = hsv2rgb(vec3(hue, 0.5, 0.25));
      
      float fresnel = pow(1.0 - abs(dot(vNormal, vec3(0,0,1))), 3.0);
      float mask = smoothstep(-0.5, 0.8, noise);
      vec3 color = rainbow * (mask * 0.6 + fresnel * 0.4);
      
      gl_FragColor = vec4(color, 1.0);
    }
  `,
  transparent: true,
  blending: THREE.AdditiveBlending,
  depthWrite: false,
  side: THREE.DoubleSide,
};

function FancyCube() {
  const spinnerRef = useRef<THREE.Group>(null);
  const ghostMaterial = useMemo(() => new THREE.ShaderMaterial(SpectralGhostShader), []);

  useFrame((state) => {
    if (spinnerRef.current) {
      spinnerRef.current.rotation.y = state.clock.elapsedTime * 0.2;
      spinnerRef.current.rotation.x = state.clock.elapsedTime * 0.1;
    }
    ghostMaterial.uniforms.uTime.value = state.clock.elapsedTime;
  });

  const size = 2.2;

  return (
    <group position={[2.5, 0, 0]}>
      <Float speed={2} rotationIntensity={0.5} floatIntensity={1}>
        <group ref={spinnerRef}>
          <mesh>
            <boxGeometry args={[size, size, size]} />
            <MeshTransmissionMaterial
              backside
              samples={16}
              resolution={1024}
              transmission={1}
              roughness={0.0}
              thickness={2.5}
              ior={1.8} // High IOR for strong refraction warping
              chromaticAberration={0.2}
              anisotropy={0.1}
              distortion={0.8}
              distortionScale={0.5}
              temporalDistortion={0.1}
              color="white"
              attenuationDistance={0.8}
              attenuationColor="#050515"
            />
            <Edges threshold={15} color="#ffffff" scale={1.01}>
               <meshBasicMaterial color="#ffffff" toneMapped={false} />
            </Edges>
          </mesh>

          <mesh scale={1.002}>
            <boxGeometry args={[size, size, size]} />
            <primitive object={ghostMaterial} attach="material" />
          </mesh>
        </group>
      </Float>
    </group>
  );
}

export default function Home() {
  // Generate random positions for star-like lightformers
  const starLightformers = useMemo(() => {
    return Array.from({ length: 50 }).map((_, i) => ({
      position: [
        (Math.random() - 0.5) * 20,
        (Math.random() - 0.5) * 20,
        (Math.random() - 0.5) * 20,
      ] as [number, number, number],
      scale: Math.random() * 0.1 + 0.05,
      intensity: Math.random() * 5 + 2,
    }));
  }, []);

  return (
    <div className="flex min-h-screen w-full flex-col md:flex-row items-center justify-center p-8 md:p-24 relative overflow-hidden bg-black">
      
      <div className="absolute inset-0 z-0">
        <Canvas dpr={[1, 2]} camera={{ position: [0, 0, 8], fov: 40 }}>
          <ambientLight intensity={0.5} />
          
          {/* BRIGHTER BACKGROUND STARS FOR REFRACTION */}
          <Stars radius={300} depth={60} count={20000} factor={12} saturation={0} fade speed={1} />
          <Sparkles scale={100} count={1000} size={4} speed={0.2} opacity={0.5} color="#ffffff" />
          
          <Suspense fallback={null}>
            <Environment resolution={1024}>
              {/* Background sphere to provide a context for the environment */}
              <mesh scale={100}>
                <sphereGeometry args={[1, 64, 64]} />
                <meshBasicMaterial color="#000000" side={THREE.BackSide} />
              </mesh>
              
              {/* STAR-LIKE LIGHTFORMERS FOR REFLECTION */}
              {starLightformers.map((star, i) => (
                <Lightformer
                  key={i}
                  form="circle"
                  intensity={star.intensity}
                  position={star.position}
                  scale={star.scale}
                  color="white"
                />
              ))}
              
              {/* Some larger lightformers for general glass sheen */}
              <Lightformer form="rect" intensity={2} rotation-x={Math.PI / 2} position={[0, 5, -5]} scale={[10, 10, 1]} />
              <Lightformer form="rect" intensity={1} rotation-y={Math.PI / 2} position={[-5, 2, 0]} scale={[10, 10, 1]} />
            </Environment>

            <FancyCube />

            <EffectComposer enableNormalPass={false}>
              <Bloom luminanceThreshold={1} mipmapBlur intensity={1.5} radius={0.4} />
            </EffectComposer>
          </Suspense>
        </Canvas>
      </div>

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

      <div className="flex-1 pointer-events-none" />
    </div>
  );
}