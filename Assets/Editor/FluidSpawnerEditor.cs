using UnityEditor;
using UnityEngine;

[CustomEditor(typeof(FluidSpawner))]
public class FluidSpawnerEditor : Editor
{
    private void OnSceneGUI()
    {
        // Extracting the component
        var spawner = (FluidSpawner)target;
        var transform = spawner.transform;

        foreach (var region in spawner.spawnRegions)
        {
            // Converting region position to world space
            var worldPos = transform.TransformPoint(region.position);
            var worldSize = region.size;

            // Only show the handle matching the current tool mode
            switch (Tools.current)
            {
                // Handler for moving the spawn region
                case Tool.Move:
                    EditorGUI.BeginChangeCheck();
                    var newWorldPos = Handles.PositionHandle(worldPos, Quaternion.identity);
                    if (EditorGUI.EndChangeCheck())
                    {
                        Undo.RecordObject(spawner, "Move Spawn Region");
                        region.position = transform.InverseTransformPoint(newWorldPos);
                    }
                    break;
                    
                // Handler for resizing the region
                case Tool.Scale:
                    EditorGUI.BeginChangeCheck();
                    var newSize = Handles.ScaleHandle(worldSize, worldPos, Quaternion.identity, 
                        HandleUtility.GetHandleSize(worldPos));
                    if (EditorGUI.EndChangeCheck())
                    {
                        Undo.RecordObject(spawner, "Resize Spawn Region");
                        region.size = newSize;
                    }
                    break;
            }
        }
    }
}