import argparse
import json
import pathlib


STANDALONE_HTML_TEMPLATE = """<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Simplicial Complex Visualizer</title>
    <style>
        body { margin: 0; overflow: hidden; background-color: #f0f0f0; font-family: sans-serif; }
        #info {
            position: absolute;
            top: 10px;
            left: 10px;
            background: rgba(255, 255, 255, 0.95);
            padding: 15px;
            border-radius: 8px;
            box-shadow: 0 2px 10px rgba(0,0,0,0.1);
            pointer-events: auto;
            width: 320px;
            max-height: 95vh;
            overflow-y: auto;
            display: flex;
            flex-direction: column;
            gap: 10px;
        }
        .control-group { border-bottom: 1px solid #eee; padding-bottom: 10px; }
        .control-group:last-child { border-bottom: none; }
        label { display: block; font-weight: bold; margin-bottom: 5px; font-size: 0.9em; color: #333; }
        textarea { width: 100%; height: 50px; padding: 5px; font-family: monospace; font-size: 0.85em; resize: vertical; box-sizing: border-box; }
        input[type="text"] { width: 65%; padding: 5px; }
        button { padding: 6px 12px; cursor: pointer; background: #007bff; color: white; border: none; border-radius: 4px; font-size: 0.9em; }
        button:hover { background: #0056b3; }
        button.secondary { background: #6c757d; }
        button.secondary:hover { background: #5a6268; }
        #stats { font-size: 0.85em; color: #666; margin-top: 5px; }
        .row { display: flex; justify-content: space-between; align-items: center; }
        #selectedCoords { background-color: #f9f9f9; color: #333; height: 80px; }
    </style>
    <script type="importmap">
        {
            "imports": {
                "three": "https://unpkg.com/three@0.160.0/build/three.module.js",
                "three/addons/": "https://unpkg.com/three@0.160.0/examples/jsm/"
            }
        }
    </script>
</head>
<body>

    <div id="info">
        <h3>Simplicial Complex</h3>
        <div id="stats">Vertices: 0 | Faces: 0</div>

        <div class="control-group">
            <label>1. Embedding (id x y z)</label>
            <textarea id="embeddingInput" placeholder="0 1.5 2.0 0.0&#10;5 0.0 1.0 1.0"></textarea>
        </div>

        <div class="control-group">
            <label>2. Faces (+/- v1 v2 v3)</label>
            <textarea id="facesInput" placeholder="+ 0 1 2&#10;- 3 4 5"></textarea>
            <button id="btnUpdate" style="width: 100%; margin-top: 5px;">Update / Render</button>
        </div>

        <div class="control-group">
            <label>3. Select Vertices (Click mesh or type)</label>
            <div class="row">
                <input type="text" id="vertexInput" placeholder="e.g. 0, 5, 12">
                <button id="btnHighlight" class="secondary">Go</button>
            </div>
        </div>

        <div class="control-group">
            <label>Selected Coordinates</label>
            <textarea id="selectedCoords" readonly placeholder="Coordinates of selected vertices will appear here..."></textarea>
        </div>
    </div>

    <script type="module">
        import * as THREE from 'three';
        import { OrbitControls } from 'three/addons/controls/OrbitControls.js';

        const INITIAL_DATA = __MESH_DATA_JSON__;

        const scene = new THREE.Scene();
        scene.background = new THREE.Color(0xf5f5f5);

        const camera = new THREE.PerspectiveCamera(45, window.innerWidth / window.innerHeight, 0.01, 1000);
        camera.position.set(0, 0, 5);
        scene.add(camera);

        const renderer = new THREE.WebGLRenderer({ antialias: true });
        renderer.setSize(window.innerWidth, window.innerHeight);
        document.body.appendChild(renderer.domElement);

        const controls = new OrbitControls(camera, renderer.domElement);
        controls.enableDamping = true;

        const ambientLight = new THREE.AmbientLight(0xffffff, 0.6);
        scene.add(ambientLight);

        const dirLight = new THREE.DirectionalLight(0xffffff, 0.6);
        dirLight.position.set(0, 0, 1);
        camera.add(dirLight);

        let currentMesh = null;
        let currentWireframe = null;
        let markersGroup = new THREE.Group();
        scene.add(markersGroup);

        let vertices = INITIAL_DATA.vertices.map(v => ({ x: v[0], y: v[1], z: v[2] }));
        let faces = INITIAL_DATA.faces.map(f => [...f]);

        const raycaster = new THREE.Raycaster();
        const mouse = new THREE.Vector2();
        let mouseDownPos = new THREE.Vector2();

        document.getElementById('btnUpdate').addEventListener('click', () => {
            applyEmbeddingUpdates();
            applyFaceUpdates();
            renderScene();
        });

        document.getElementById('btnHighlight').addEventListener('click', () => {
            const input = document.getElementById('vertexInput').value;
            highlightVertices(input);
        });

        window.addEventListener('mousedown', (event) => {
            if (event.target !== renderer.domElement) return;
            mouseDownPos.set(event.clientX, event.clientY);
        });

        window.addEventListener('mouseup', (event) => {
            if (event.target !== renderer.domElement) return;
            const dist = mouseDownPos.distanceTo(new THREE.Vector2(event.clientX, event.clientY));
            if (dist < 5) {
                onMouseClick(event);
            }
        });

        function onMouseClick(event) {
            if (!currentMesh) return;

            mouse.x = (event.clientX / window.innerWidth) * 2 - 1;
            mouse.y = -(event.clientY / window.innerHeight) * 2 + 1;

            raycaster.setFromCamera(mouse, camera);
            const intersects = raycaster.intersectObject(currentMesh);

            if (intersects.length > 0) {
                const intersection = intersects[0];
                const face = intersection.face;
                const point = intersection.point;
                const indices = [face.a, face.b, face.c];

                const posAttribute = currentMesh.geometry.attributes.position;
                let minDst = Infinity;
                let closestVertexId = -1;

                indices.forEach(idx => {
                    const vPos = new THREE.Vector3(
                        posAttribute.getX(idx),
                        posAttribute.getY(idx),
                        posAttribute.getZ(idx)
                    );
                    const dst = point.distanceTo(vPos);
                    if (dst < minDst) {
                        minDst = dst;
                        closestVertexId = idx;
                    }
                });

                if (closestVertexId !== -1) {
                    toggleVertexSelection(closestVertexId);
                }
            }
        }

        function applyEmbeddingUpdates() {
            const text = document.getElementById('embeddingInput').value.trim();
            if (!text) return;

            const lines = text.split('\\n');
            lines.forEach(line => {
                const parts = line.trim().split(/\\s+/);
                if (parts.length < 4) return;
                const id = parseInt(parts[0]);
                const x = parseFloat(parts[1]);
                const y = parseFloat(parts[2]);
                const z = parseFloat(parts[3]);
                if (isNaN(id) || isNaN(x)) return;

                if (id >= vertices.length) {
                    for (let k = vertices.length; k <= id; k++) {
                        vertices[k] = { x: 0, y: 0, z: 0 };
                    }
                }
                vertices[id] = { x, y, z };
            });
        }

        function applyFaceUpdates() {
            const text = document.getElementById('facesInput').value.trim();
            if (!text) return;

            const lines = text.split('\\n');
            lines.forEach(line => {
                const parts = line.trim().split(/\\s+/);
                if (parts.length < 4) return;
                const op = parts[0];
                const idxs = [parseInt(parts[1]), parseInt(parts[2]), parseInt(parts[3])];
                const sortedIdxs = [...idxs].sort((a, b) => a - b);

                const existingIndex = faces.findIndex(f => {
                    const sortedF = [...f].sort((a, b) => a - b);
                    return sortedF[0] === sortedIdxs[0] && sortedF[1] === sortedIdxs[1] && sortedF[2] === sortedIdxs[2];
                });

                if (op === '+' && existingIndex === -1) {
                    faces.push(idxs);
                } else if (op === '-' && existingIndex !== -1) {
                    faces.splice(existingIndex, 1);
                }
            });
        }

        function renderScene() {
            document.getElementById('stats').innerText = `Vertices: ${vertices.length} | Faces: ${faces.length}`;

            if (currentMesh) {
                scene.remove(currentMesh);
                currentMesh.geometry.dispose();
                currentMesh.material.dispose();
            }
            if (currentWireframe) {
                scene.remove(currentWireframe);
                currentWireframe.geometry.dispose();
                currentWireframe.material.dispose();
            }

            clearMarkers();
            document.getElementById('vertexInput').value = "";
            document.getElementById('selectedCoords').value = "";

            if (vertices.length === 0) return;

            const positions = [];
            vertices.forEach(v => positions.push(v.x, v.y, v.z));

            const indices = [];
            faces.forEach(f => {
                if (f[0] < vertices.length && f[1] < vertices.length && f[2] < vertices.length) {
                    indices.push(f[0], f[1], f[2]);
                }
            });

            const geometry = new THREE.BufferGeometry();
            geometry.setAttribute('position', new THREE.Float32BufferAttribute(positions, 3));
            geometry.setIndex(indices);
            geometry.center();

            const material = new THREE.MeshPhongMaterial({
                color: 0xE0F7FA,
                side: THREE.DoubleSide,
                flatShading: true,
                shininess: 0,
                specular: 0x000000
            });

            currentMesh = new THREE.Mesh(geometry, material);
            scene.add(currentMesh);

            const wireGeo = new THREE.WireframeGeometry(geometry);
            const wireMat = new THREE.LineBasicMaterial({ color: 0x000000, opacity: 0.25, transparent: true });
            currentWireframe = new THREE.LineSegments(wireGeo, wireMat);
            scene.add(currentWireframe);

            geometry.computeBoundingSphere();
            const r = geometry.boundingSphere.radius || 1;
            camera.position.set(0, 0, r * 2.5);
            controls.target.set(0, 0, 0);
            controls.update();
        }

        function toggleVertexSelection(id) {
            const inputField = document.getElementById('vertexInput');
            let currentIds = inputField.value.split(/[\\s,]+/)
                .map(s => parseInt(s))
                .filter(n => !isNaN(n));

            const idx = currentIds.indexOf(id);
            if (idx === -1) {
                currentIds.push(id);
            } else {
                currentIds.splice(idx, 1);
            }

            const newValue = currentIds.join(', ');
            inputField.value = newValue;
            highlightVertices(newValue);
        }

        function highlightVertices(inputString) {
            clearMarkers();
            if (!currentMesh) return;

            const posAttr = currentMesh.geometry.attributes.position;
            const count = posAttr.count;

            const ids = inputString.split(/[\\s,]+/)
                .map(s => parseInt(s))
                .filter(n => !isNaN(n) && n >= 0 && n < count);

            currentMesh.geometry.computeBoundingSphere();
            const scale = (currentMesh.geometry.boundingSphere.radius || 1) * 0.005;

            const coordOutput = [];

            ids.forEach(id => {
                const x = posAttr.getX(id);
                const y = posAttr.getY(id);
                const z = posAttr.getZ(id);

                const marker = new THREE.Mesh(
                    new THREE.SphereGeometry(scale, 16, 16),
                    new THREE.MeshBasicMaterial({ color: 0xff0000 })
                );
                marker.position.set(x, y, z);
                markersGroup.add(marker);

                const sprite = makeTextSprite(id.toString(), { fontsize: 64 });
                sprite.position.set(x, y + scale * 8.0, z);
                sprite.scale.set(scale * 25, scale * 12.5, 1);
                markersGroup.add(sprite);

                if (vertices[id]) {
                    const v = vertices[id];
                    coordOutput.push(`${id} ${v.x} ${v.y} ${v.z}`);
                }
            });

            document.getElementById('selectedCoords').value = coordOutput.join('\\n');
        }

        function clearMarkers() {
            while (markersGroup.children.length > 0) {
                const child = markersGroup.children[0];
                markersGroup.remove(child);
                if (child.geometry) child.geometry.dispose();
                if (child.material) {
                    if (child.material.map) child.material.map.dispose();
                    child.material.dispose();
                }
            }
        }

        function makeTextSprite(message, parameters) {
            if (parameters === undefined) parameters = {};
            const fontface = parameters.hasOwnProperty("fontface") ? parameters["fontface"] : "Arial";
            const fontsize = parameters.hasOwnProperty("fontsize") ? parameters["fontsize"] : 48;

            const canvas = document.createElement('canvas');
            const context = canvas.getContext('2d');
            canvas.width = 512;
            canvas.height = 256;

            context.font = fontsize + "px " + fontface;
            context.fillStyle = "rgba(0, 0, 0, 1.0)";
            context.textAlign = "center";
            context.textBaseline = "middle";
            context.fillText(message, canvas.width / 2, canvas.height / 2);

            const texture = new THREE.Texture(canvas);
            texture.needsUpdate = true;
            const spriteMaterial = new THREE.SpriteMaterial({ map: texture });
            return new THREE.Sprite(spriteMaterial);
        }

        function animate() {
            requestAnimationFrame(animate);
            controls.update();
            renderer.render(scene, camera);
        }

        window.addEventListener('resize', () => {
            camera.aspect = window.innerWidth / window.innerHeight;
            camera.updateProjectionMatrix();
            renderer.setSize(window.innerWidth, window.innerHeight);
        });

        renderScene();
        animate();
    </script>
</body>
</html>
"""


def parse_competitor_file(filepath):
    """Parses a competitor .txt file into a JSON-serializable dict with vertices and triangular faces."""
    with open(filepath, "r") as f:
        tokens = f.read().split()

    if len(tokens) < 2:
        raise ValueError(f"Invalid competitor file: {filepath}")

    ptr = 0
    num_vertices = int(tokens[ptr])
    ptr += 1
    num_faces = int(tokens[ptr])
    ptr += 1

    vertices = []
    for _ in range(num_vertices):
        x = float(tokens[ptr])
        y = float(tokens[ptr + 1])
        z = float(tokens[ptr + 2])
        ptr += 3
        vertices.append([x, y, z])

    faces = []
    for _ in range(num_faces):
        k = int(tokens[ptr])
        ptr += 1
        if k == 3:
            faces.append([int(tokens[ptr]), int(tokens[ptr + 1]), int(tokens[ptr + 2])])
        ptr += k

    return {"vertices": vertices, "faces": faces}


def write_competitor_html(mesh_data, output_file):
    """Generates a standalone interactive HTML file for the given competitor mesh."""
    html_content = STANDALONE_HTML_TEMPLATE.replace("__MESH_DATA_JSON__", json.dumps(mesh_data))
    output_path = pathlib.Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w") as f:
        f.write(html_content)


def main():
    parser = argparse.ArgumentParser(description="Visualize a competitor simplicial complex.")
    parser.add_argument("--in", dest="input_file", required=True, help="Input competitor .txt file.")
    parser.add_argument("--out", dest="output_file", required=True, help="Output standalone .html file.")
    args = parser.parse_args()

    mesh_data = parse_competitor_file(args.input_file)
    write_competitor_html(mesh_data, args.output_file)


if __name__ == "__main__":
    main()
