import { StrictMode } from 'react'
import { createRoot } from 'react-dom/client'
import './index.css'
import App from './App.tsx'
import Viewer from './Viewer.tsx'

// The server serves this same bundle at every non-API path (backend/app/main.py's
// mount_frontend); /view is the logged-out page, everything else the full app (#372).
const isViewer = window.location.pathname.replace(/\/+$/, '') === '/view'

createRoot(document.getElementById('root')!).render(
  <StrictMode>
    {isViewer ? <Viewer /> : <App />}
  </StrictMode>,
)
