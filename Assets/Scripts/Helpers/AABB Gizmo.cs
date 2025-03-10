using UnityEngine;

public class AABBGizmo : MonoBehaviour
{
    [Header("Gizmo")] 
    [SerializeField] public Color gizmoColor;
    
    private void OnDrawGizmos()
    {
        // Drawing the simulation bounding box
        Gizmos.color = gizmoColor;
        Gizmos.DrawWireCube(transform.position, transform.localScale);
    }
}
