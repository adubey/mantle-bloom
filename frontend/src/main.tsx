import { StrictMode } from 'react'
import { createRoot } from 'react-dom/client'
import './index.css'
import App from './App.tsx'
import Viewer from './Viewer.tsx'

// The server serves this same bundle at every non-API path (backend/app/main.py's
// mount_frontend); /view is the logged-out page, with /demo/view for a local preview.
const viewerPath = window.location.pathname.replace(/\/+$/, '')
const isViewer = viewerPath === '/view' || viewerPath === '/demo/view'

createRoot(document.getElementById('root')!).render(
  <StrictMode>
    {isViewer ? <Viewer /> : <App />}
  </StrictMode>,
)
