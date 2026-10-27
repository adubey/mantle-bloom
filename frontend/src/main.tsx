import { StrictMode } from 'react'
import { createRoot } from 'react-dom/client'
import './index.css'
import App from './App.tsx'
import Viewer from './Viewer.tsx'
import NotFound from './NotFound.tsx'

// The server serves this same bundle at every non-API path (backend/app/main.py's
// mount_frontend); /view is the logged-out page, with /demo/view for a local preview.
const path = window.location.pathname.replace(/\/+$/, '') || '/'
const isViewer = path === '/view' || path === '/demo/view'
const isApp = path === '/'

createRoot(document.getElementById('root')!).render(
  <StrictMode>
    {isViewer ? <Viewer /> : isApp ? <App /> : <NotFound />}
  </StrictMode>,
)
